#include "vertex_data.h"

#include <Eigen/Dense>
#include <fmt/format.h>
#include <iostream>
#include <solp.h>
#include <numeric>

template<int d>
static Eigen::Matrix<double, d * d, d * d> kronecker_prod(const Eigen::Matrix<double, d, d> &a,
                                                          const Eigen::Matrix<double, d, d> &b) {
	Eigen::Matrix<double, d * d, d * d> res;
	for(int i1 = 0; i1 < d; i1++) {
		for(int j1 = 0; j1 < d; j1++) {
			for(int i2 = 0; i2 < d; i2++) {
				for(int j2 = 0; j2 < d; j2++) {
					int i = d * i1 + i2;
					int j = d * j1 + j2;
					res(i, j) = a(i1, j1) * b(i2, j2);
				}
			}
		}
	}

	return res;
}

using mat = Eigen::Matrix<double, jm::basis_size, jm::basis_size>;

static auto spinops(int nspinhalfs) {
	std::array<mat, 3> S;
	S[0] = mat::Zero();
	S[1] = mat::Zero();
	S[2] = mat::Zero();
	
	if(nspinhalfs == 1) {
		S[0].diagonal(1) << 1, 0, 0;
		S[2].diagonal() << 0.5, -0.5, 0, 0;
	} else if(nspinhalfs == 2) {
		S[0].diagonal(1) << 0, sqrt(2), sqrt(2);
		S[2].diagonal() << 0, 1, 0, -1;
	} else {
		assert(false);
	}

	S[1] = S[0].transpose();

	return S;
}

static mat jin_term(int nspinhalfs) {
	auto S = spinops(nspinhalfs);

	mat S2 = 0.5 * (S[0] * S[1] + S[1] * S[0]) + S[2] * S[2];

	if(nspinhalfs == 1) {
		return mat::Zero();
	} else if(nspinhalfs == 2) {
		return 0.5*(S2 - mat::Identity()*1.5);
	}
	assert(false);
	return mat::Zero();
}

void vertex_data::construct_vertices(const bond &b, const site &si, const site &sj) {
	std::vector<mat> S(3);

	auto Si = spinops(si.nspinhalfs);
	auto Sj = spinops(sj.nspinhalfs);

	Eigen::Matrix<double, jm::basis_size * jm::basis_size, jm::basis_size *jm::basis_size> H =
	    b.J * (0.5 * (kronecker_prod(Si[0], Sj[1]) + kronecker_prod(Si[1], Sj[0])) +
	           kronecker_prod(Si[2], Sj[2]));
	mat Id = mat::Identity();
	H += si.Jin * kronecker_prod(jin_term(si.nspinhalfs),Id);
	H += sj.Jin * kronecker_prod(Id,jin_term(sj.nspinhalfs));

	H.diagonal() *= -1; // from -β ;)

	double epsilon = 0.;
	energy_offset = H.diagonal().minCoeff()-epsilon;
	H -= decltype(H)::Identity()*energy_offset;

	for(uint8_t l0 = 0; l0 < jm::basis_size; l0++) {
		for(uint8_t l1 = 0; l1 < jm::basis_size; l1++) {
			for(uint8_t l2 = 0; l2 < jm::basis_size; l2++) {
				for(uint8_t l3 = 0; l3 < jm::basis_size; l3++) {
					double w = H(l0 * jm::basis_size + l1, l2 * jm::basis_size + l3);
					int lbi = jm::local_basis_size(si);
					int lbj = jm::local_basis_size(sj);
					if(w != 0 && l0 < lbi && l1 < lbj && l2 < lbi && l3 < lbj) {
						legstates.push_back({jm{l0}, jm{l1}, jm{l2}, jm{l3}});
						weights_.push_back(w);
						assert(w >= 0);
					}
				}
			}
		}
	}
}

int vertex_data::vertex_change_apply(int vertex, int leg_in, jm_action action_in, int leg_out,
                                     jm_action action_out) const {
	auto legstate = legstates.at(vertex);
	legstate[leg_in].apply(action_in);
	legstate[leg_out].apply(action_out);

	auto it = std::find(legstates.begin(), legstates.end(), legstate);
	if(it == legstates.end()) {
		return -1;
	}

	return it - legstates.begin();
}

void vertex_data::init_code_to_idx() {
	for(size_t i = 0; i < legstates.size(); i++) {
		const auto &ls = legstates[i];
		code_to_idx_[opercode::make_vertex(0,ls[0],ls[1],ls[2],ls[3]).vertex()] = i;
	}
}
	

vertex_data::vertex_data(const bond &b, const site &si, const site &sj) {
	construct_vertices(b, si, sj);
	transitions_.resize(legstates.size() * jm::basis_size * leg_count);

	struct vertex_change {
		int leg_in{};
		jm_action action_in{};
		int leg_out{};
		jm_action action_out{};

		vertex_change inverse() const {
			return vertex_change{leg_out, action_out, leg_in, action_in};
		}
		bool operator==(const vertex_change &other) const {
			return leg_in == other.leg_in && action_in == other.action_in &&
			       leg_out == other.leg_out && action_out == other.action_out;
		}
	};

	std::vector<vertex_change> variables;
	for(jm_action action_in = 0; action_in < jm::basis_size; action_in++) {
		for(int leg_in = 0; leg_in < leg_count; leg_in++) {
			for(jm_action action_out = 0; action_out < jm::basis_size; action_out++) {
				for(int leg_out = 0; leg_out < leg_count; leg_out++) {
					vertex_change vc{leg_in, action_in, leg_out, action_out};

					auto it = std::find_if(variables.begin(), variables.end(),
					                       [&](auto x) { return vc.inverse() == x; });
					if(it == variables.end()) {
						variables.push_back(vc);
					}
				}
			}
		}
	}

	std::vector<solp::constraint> constraints;
	std::vector<double> objective(variables.size(), 0);
	std::transform(
	    variables.begin(), variables.end(), objective.begin(), [](const vertex_change &vc) {
		    bool exclude_bounce = vc.leg_in == vc.leg_out && vc.action_in == vc.action_out;
		    bool exclude_noop = vc.action_in == 0 && vc.action_out == 0;
		    return exclude_bounce || exclude_noop;
	    });

	for(jm_action action_out = 0; action_out < jm::basis_size; action_out++) {
		for(int leg_out = 0; leg_out < leg_count; leg_out++) {
			std::vector<double> coeff(variables.size(), 0);

			for(size_t i = 0; i < variables.size(); i++) {
				coeff[i] =
				    (variables[i].leg_in == leg_out && variables[i].action_in == action_out) ||
				    (variables[i].leg_out == leg_out && variables[i].action_out == action_out);
			}
			constraints.push_back(solp::constraint{coeff, 0});
		}
	}

	for(size_t v = 0; v < weights_.size(); v++) {
		for(jm_action action_in = 0; action_in < jm::basis_size; action_in++) {
			for(int leg_in = 0; leg_in < leg_count; leg_in++) {
				std::vector<double> ws(leg_count * jm::basis_size); // [action_out*leg_count+leg_out]
				std::vector<int> targets(leg_count * jm::basis_size);

				for(jm_action action_out = 0; action_out < jm::basis_size; action_out++) {
					for(int leg_out = 0; leg_out < leg_count; leg_out++) {
						int target = vertex_change_apply(v, leg_in, action_in, leg_out, action_out);
						int out = action_out * leg_count + leg_out;

						targets[out] = target;
						constraints[out].rhs = target >= 0 ? weights_[target] : 0;
					}
				}

				solp::result result = solp::solve(objective, constraints, solp::options{1e-10});

				for(size_t i = 0; i < variables.size(); i++) {
					const auto &var = variables[i];
					int in = var.action_in * leg_count + var.leg_in;
					int out = var.action_out * leg_count + var.leg_out;

					if(targets[in] >= 0) {
						transitions_[targets[in] * jm::basis_size * leg_count + in].targets[out] =
						    targets[out];
						double norm = constraints[in].rhs == 0 ? 1 : constraints[in].rhs;
						transitions_[targets[in] * jm::basis_size * leg_count + in].probs[out] =
						    result.x[i] > 0 ? result.x[i]/norm : 0;
					}

					in = var.inverse().action_in * leg_count + var.inverse().leg_in;
					out = var.inverse().action_out * leg_count + var.inverse().leg_out;

					if(targets[in] >= 0) {
						transitions_[targets[in] * jm::basis_size * leg_count + in].targets[out] =
						    targets[out];
						double norm = constraints[in].rhs == 0 ? 1 : constraints[in].rhs;
						transitions_[targets[in] * jm::basis_size * leg_count + in].probs[out] =
						    result.x[i] > 0 ? result.x[i]/norm : 0;
					}
				}
			}
		}
	}

	for(auto &t : transitions_) {
		std::partial_sum(t.probs.begin(), t.probs.end(), t.probs.begin());
		assert(!t.invalid());
	}

	init_code_to_idx();
}

void vertex_data::print(const site &si, const site &sj) const {
	for(size_t v = 0; v < legstates.size(); v++) {
		for(size_t in = 0; in < jm::basis_size * leg_count; in++) {
			const auto &trans = transitions_[v * jm::basis_size * leg_count + in];
			std::cout << fmt::format("{:.2f}|{:2d}\n",
			                         fmt::join(trans.probs.begin(), trans.probs.end(), " "),
			                         fmt::join(trans.targets.begin(), trans.targets.end(), " "));
		}
		std::cout << "\n";
	}

	std::cout << "\nlegstates:\n";

	int idx{};
	for(const auto &ls : legstates) {
		int idxtest = code_to_idx_.at(opercode::make_vertex(0,ls[0],ls[1],ls[2],ls[3]).vertex());
		std::cout << fmt::format("{}({}): [{}{} {}{}] ~ {}\n", idx, idxtest, ls[0].name(si), ls[1].name(sj), ls[2].name(si), ls[3].name(sj), weights_[idx]);
		idx++;
	}

}

