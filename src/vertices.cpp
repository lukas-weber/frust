#include "vertices.h"

#include <Eigen/Dense>
#include <fmt/format.h>
#include <iostream>
#include <solp.h>

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

static auto spinops(int nspinhalfs) {
	std::array<Eigen::Matrix<double, vertex_data::basis_size, vertex_data::basis_size>, 3> S;
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

void vertex_data::construct_vertices(const bond &b, const site &si, const site &sj) {
	using mat = Eigen::Matrix<double, basis_size, basis_size>;

	std::vector<mat> S(3);

	auto Si = spinops(si.nspinhalfs);
	auto Sj = spinops(sj.nspinhalfs);

	Eigen::Matrix<double, basis_size * basis_size, basis_size *basis_size> H =
	    b.J * (0.5 * (kronecker_prod(Si[0], Sj[1]) + kronecker_prod(Si[1], Sj[0])) +
	           kronecker_prod(Si[2], Sj[2]));
	mat S2i = 0.5 * (Si[0] * Si[1] + Si[1] * Si[0]) + Si[2] * Si[2];
	mat S2j = 0.5 * (Sj[0] * Sj[1] + Sj[1] * Sj[0]) + Sj[2] * Sj[2];
	mat Id = mat::Identity();
	H += b.Ki * kronecker_prod(S2i, Id);
	H += b.Kj * kronecker_prod(Id, S2j);

	for(uint32_t l0 = 0; l0 < basis_size; l0++) {
		for(uint32_t l1 = 0; l1 < basis_size; l1++) {
			for(uint32_t l2 = 0; l2 < basis_size; l2++) {
				for(uint32_t l3 = 0; l3 < basis_size; l3++) {
					double w = H(l0 * basis_size + l1, l2 * basis_size + l3);
					if(w != 0) {
						legstates_.push_back({jm{l0}, jm{l1}, jm{l2}, jm{l3}});
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
	auto legstate = legstates_.at(vertex);
	legstate[leg_in].apply(action_in);
	legstate[leg_out].apply(action_out);

	auto it = std::find(legstates_.begin(), legstates_.end(), legstate);
	if(it == legstates_.end()) {
		return -1;
	}

	return it - legstates_.begin();
}

vertex_data::vertex_data(const bond &b, const site &si, const site &sj) {
	construct_vertices(b, si, sj);
	transitions.resize(legstates_.size() * basis_size * leg_count);

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
	for(jm_action action_in = 0; action_in < basis_size; action_in++) {
		for(int leg_in = 0; leg_in < leg_count; leg_in++) {
			for(jm_action action_out = 0; action_out < basis_size; action_out++) {
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

	for(jm_action action_out = 0; action_out < basis_size; action_out++) {
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
		std::cout << fmt::format("v = {}/{}\n", v, weights_.size());
		for(jm_action action_in = 0; action_in < basis_size; action_in++) {
			for(int leg_in = 0; leg_in < leg_count; leg_in++) {
				std::vector<double> ws(leg_count * basis_size); // [action_out*leg_count+leg_out]
				std::vector<int> targets(leg_count * basis_size);

				for(jm_action action_out = 0; action_out < basis_size; action_out++) {
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
						transitions[targets[in] * basis_size * leg_count + in].targets[out] =
						    targets[out];
						double norm = constraints[in].rhs == 0 ? 1 : constraints[in].rhs;
						transitions[targets[in] * basis_size * leg_count + in].probs[out] =
						    result.x[i] > 0 ? result.x[i]/norm : 0;
					}

					in = var.inverse().action_in * leg_count + var.inverse().leg_in;
					out = var.inverse().action_out * leg_count + var.inverse().leg_out;

					if(targets[in] >= 0) {
						transitions[targets[in] * basis_size * leg_count + in].targets[out] =
						    targets[out];
						double norm = constraints[in].rhs == 0 ? 1 : constraints[in].rhs;
						transitions[targets[in] * basis_size * leg_count + in].probs[out] =
						    result.x[i] > 0 ? result.x[i]/norm : 0;
					}
				}
			}
		}
	}
}

void vertex_data::print() const {
	for(size_t v = 0; v < legstates_.size(); v++) {
		for(size_t in = 0; in < basis_size * leg_count; in++) {
			const auto &trans = transitions[v * basis_size * leg_count + in];
			std::cout << fmt::format("{:.2f}|{:2d}\n",
			                         fmt::join(trans.probs.begin(), trans.probs.end(), " "),
			                         fmt::join(trans.targets.begin(), trans.targets.end(), " "));
		}
		std::cout << "\n";
	}
}
