#include "vertex_data.h"

#include <Eigen/Dense>
#include <fmt/format.h>
#include <iostream>
#include <solp.h>
#include <numeric>
#include "util.h"

static Eigen::MatrixXd onsite_term(const uc_site &s) {
	const auto &b = s.basis;

	Eigen::MatrixXd res = res.Zero(b.size(), b.size());

	for(int i = 0; i < b.nspinhalfs; i++) {
		auto spini = b.spinop(i);
		for(int j = 0; j < i; j++) {
			auto spinj = b.spinop(j);
			res += 0.5*(spini[0]*spinj[1] + spini[1]*spinj[0]) + spini[2]*spinj[2];
		}
	}
	return s.Jin * res;
}

static Eigen::MatrixXd bond_term(const uc_bond &b, const uc_site &si, const uc_site &sj) {
	assert(static_cast<int>(b.J.size()) == si.basis.nspinhalfs*sj.basis.nspinhalfs);

	int dim = si.basis.size()*sj.basis.size();
	Eigen::MatrixXd res = Eigen::MatrixXd::Zero(dim, dim);
	
	for(int i = 0; i < si.basis.nspinhalfs; i++) {
		auto spini = si.basis.spinop(i);
		for(int j = 0; j < sj.basis.nspinhalfs; j++) {
			auto spinj = sj.basis.spinop(j);
			res += b.J[sj.basis.nspinhalfs*i + j] * scalar_product(spini, spinj);
		}
	}

	return res;
}

void vertex_data::construct_vertices(const uc_bond &b, const uc_site &si, const uc_site &sj, double tolerance) {
	int dim = si.basis.size()*sj.basis.size();
	Eigen::MatrixXd H = Eigen::MatrixXd::Zero(dim, dim);
	auto Idi = Eigen::MatrixXd::Identity(si.basis.size(),si.basis.size());
	auto Idj = Eigen::MatrixXd::Identity(sj.basis.size(),sj.basis.size());
	H += kronecker_prod(onsite_term(si),Idj)/si.coordination;
	H += kronecker_prod(Idi,onsite_term(sj))/sj.coordination;
	H += bond_term(b, si, sj);

	H *= -1; // from -β ;)


	double epsilon = fabs(H.maxCoeff())/2;
	energy_offset = H.diagonal().minCoeff()-epsilon;
	H -= decltype(H)::Identity(dim, dim)*energy_offset;

	assert(H.isApprox(H.transpose()));

	for(state_idx l0 = 0; l0 < si.basis.size(); l0++) {
		for(state_idx l1 = 0; l1 < sj.basis.size(); l1++) {
			for(state_idx l2 = 0; l2 < si.basis.size(); l2++) {
				for(state_idx l3 = 0; l3 < sj.basis.size(); l3++) {
					double w = H(l0 * sj.basis.size() + l1, l2 * sj.basis.size() + l3);
					if(fabs(w) > tolerance) {
						legstates_.push_back({l0, l1, l2, l3});
						weights_.push_back(fabs(w));
						signs_.push_back(w >= 0 ? 1 : -1);

						if(l0 == l2 && l1 == l3) {
							diagonal_vertices_[l0*site_basis::max_size + l1] = vertexcode{true, static_cast<uint32_t>(weights_.size()-1)};
						}
					}
				}
			}
		}
	}
}

vertexcode vertex_data::wrap_vertex_idx(int vertex_idx) {
	if(vertex_idx < 0) {
		return vertexcode{};
	}
	
	const auto &ls = legstates_[vertex_idx];
	return vertexcode{ls[0]==ls[2] && ls[1] == ls[3], static_cast<uint32_t>(vertex_idx)};
}

int vertex_data::vertex_change_apply(const site_basis &bi, const site_basis &bj, int vertex, int leg_in, worm_idx worm_in_idx, int leg_out, worm_idx worm_out_idx) const {
	auto legstate = legstates_.at(vertex);
	
	const auto &basis_in = leg_in&1 ? bj : bi;
	const auto &basis_out = leg_out&1 ? bj : bi;

	const auto &worm_in =  basis_in.worms[worm_in_idx];
	const auto &worm_out = basis_out.worms[worm_out_idx];
	

	legstate[leg_in] = worm_in.action[legstate[leg_in]];
	if(legstate[leg_in] == site_basis::invalid) {
		if(leg_in == leg_out && worm_in_idx == worm_out.inverse_idx) {
			return vertex;
		}
		return -1;
	}
		
	legstate[leg_out] = worm_out.action[legstate[leg_out]];

	auto it = std::find(legstates_.begin(), legstates_.end(), legstate);
	if(it == legstates_.end()) {
		return -1;
	}

	return it - legstates_.begin();
}

vertex_data::vertex_data(const uc_bond &b, const uc_site &si, const uc_site &sj) {
	const double tolerance = 1e-10;
	construct_vertices(b, si, sj, tolerance);
	transitions_.resize(legstates_.size() * site_basis::worm_count * leg_count);

	struct vertex_change {
		const site_basis &bi;
		const site_basis &bj;

		int leg_in{};
		worm_idx worm_in{};
		int leg_out{};
		worm_idx worm_out{};

		vertex_change inverse() const {
			const auto &basis_in = leg_in&1 ? bj : bi;
			const auto &basis_out = leg_out&1 ? bj : bi;
			return vertex_change{bi, bj, leg_out, basis_out.worms[worm_out].inverse_idx, leg_in, basis_in.worms[worm_in].inverse_idx};
		}
		bool operator==(const vertex_change &other) const {
			return leg_in == other.leg_in && worm_in == other.worm_in &&
			       leg_out == other.leg_out && worm_out == other.worm_out;
		}
	};

	std::vector<vertex_change> variables;
	for(worm_idx worm_in = 0; worm_in < site_basis::worm_count; worm_in++) {
		for(int leg_in = 0; leg_in < leg_count; leg_in++) {
			for(worm_idx worm_out = 0; worm_out < site_basis::worm_count; worm_out++) {
				for(int leg_out = 0; leg_out < leg_count; leg_out++) {
					vertex_change vc{si.basis, sj.basis, leg_in, worm_in, leg_out, worm_out};

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
	    variables.begin(), variables.end(), objective.begin(), [&](const vertex_change &vc) {
		const auto &site_out = vc.leg_out&1 ? sj : si;
		bool exclude_bounce = vc.leg_in == vc.leg_out && vc.worm_in == site_out.basis.worms[vc.worm_out].inverse_idx;
		return exclude_bounce;
	    });

	for(int worm_out = 0; worm_out < site_basis::worm_count; worm_out++) {
		for(int leg_out = 0; leg_out < leg_count; leg_out++) {
			std::vector<double> coeff(variables.size(), 0);

			for(size_t i = 0; i < variables.size(); i++) {
				coeff[i] =
				    (variables[i].leg_in == leg_out && variables[i].worm_in == worm_out) ||
				    (variables[i].inverse().leg_in == leg_out && variables[i].inverse().worm_in == worm_out);
			}
			constraints.push_back(solp::constraint{coeff, 0});
		}
	}

	for(auto &t : transitions_) {
		t.probs.fill(-1);
	}

	for(size_t v = 0; v < weights_.size(); v++) {
		for(int func_in = 0; func_in < site_basis::worm_count; func_in++) {
			for(int leg_in = 0; leg_in < leg_count; leg_in++) {
				std::vector<double> ws(leg_count * site_basis::worm_count); // [func_out*leg_count+leg_out]
				std::vector<int> targets(leg_count * site_basis::worm_count);

				for(int func_out = 0; func_out < site_basis::worm_count; func_out++) {
					for(int leg_out = 0; leg_out < leg_count; leg_out++) {
						int target = vertex_change_apply(si.basis, sj.basis, v, leg_in, func_in, leg_out, func_out);
						int out = func_out * leg_count + leg_out;

						targets[out] = target;
						constraints[out].rhs = target >= 0 ? weights_[target] : 0;
					}
				}
				solp::result result = solp::solve(objective, constraints, solp::options{tolerance});

				for(size_t i = 0; i < variables.size(); i++) {
					const auto &var = variables[i];
					int in = var.worm_in * leg_count + var.leg_in;
					int out = var.worm_out * leg_count + var.leg_out;
					int in_inv = var.inverse().worm_in * leg_count + var.inverse().leg_in;

					assert(result.x[i] >= -tolerance);
					if(result.x[i] < tolerance) {
						result.x[i] = 0;
					}

					if(targets[in] >= 0) {
						transitions_[targets[in] * site_basis::worm_count * leg_count + in].targets[out] =
						    wrap_vertex_idx(targets[in_inv]);
						double norm = constraints[in].rhs == 0 ? 1 : constraints[in].rhs;
						assert(norm > 0);
						transitions_[targets[in] * site_basis::worm_count * leg_count + in].probs[out] =
						    result.x[i]/norm;
					}

					in = var.inverse().worm_in * leg_count + var.inverse().leg_in;
					out = var.inverse().worm_out * leg_count + var.inverse().leg_out;
					in_inv = var.worm_in * leg_count + var.leg_in;

					if(targets[in] >= 0) {
						transitions_[targets[in] * site_basis::worm_count * leg_count + in].targets[out] =
						    wrap_vertex_idx(targets[in_inv]);
						double norm = constraints[in].rhs == 0 ? 1 : constraints[in].rhs;
						assert(norm > 0);
						transitions_[targets[in] * site_basis::worm_count * leg_count + in].probs[out] =
						    result.x[i]/norm;
					}
				}
			}
		}
	}

	for(auto &t : transitions_) {
		int idx{};
		for(auto &p : t.probs) {
			if(p > 0 && t.targets[idx].invalid()) {
				throw std::runtime_error{fmt::format("it is possible to reach an invalid vertex with p={}", p)};
			}
			idx++;
		}

		std::partial_sum(t.probs.begin(), t.probs.end(), t.probs.begin());
		assert(t.probs.back() < 1.0000000001);
		assert(!t.invalid());
	}
}

void vertex_data::transition::print() const {
		auto tmp = probs;
		for(size_t i = 1; i < probs.size(); i++) {
			tmp[i] = probs[i]-probs[i-1];
			if(tmp[i] > -1e-8) {
				tmp[i] = fabs(tmp[i]);
			}
		}
		//std::cout << fmt::format("{:.2f}|{:2d}\n",
		//                         fmt::join(tmp.begin(), tmp.end(), " "),
		//                         fmt::join(targets.begin(), targets.end(), " "));
}
	

void vertex_data::print(const site_basis &bi, const site_basis &bj) const {
	for(size_t v = 0; v < legstates_.size(); v++) {
		for(size_t in = 0; in < site_basis::worm_count * leg_count; in++) {
			const auto &trans = transitions_[v * site_basis::worm_count * leg_count + in];
			trans.print();
		}
		std::cout << "\n";
	}

	std::cout << "\nlegstates:\n";
	std::vector<int> idxs(weights_.size());
	int idx{};
	for(const auto &ls : legstates_) {
		(void)(ls);
		idxs[idx] = idx;
		idx++;
	}
	//std::sort(idxs.begin(), idxs.end(), [&](int i, int j) {return weights_[i] < weights_[j];});


	for(const auto &idx : idxs) {
		const auto &ls = legstates_[idx];
		std::cout << fmt::format("{}({}): [{}{} {}{}] ~ {}\n", idx, signs_[idx] > 0 ? '+' : '-', bi.states[ls[0]].name, bj.states[ls[1]].name, bi.states[ls[2]].name, bj.states[ls[3]].name, weights_[idx]);
	}

}

