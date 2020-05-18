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

	if(static_cast<int>(s.Jin.size()) != b.nspinhalfs*(b.nspinhalfs-1)/2) {
		throw std::runtime_error{"Jin does not have the right shape"};
	}

	int idx = 0;
	for(int i = 0; i < b.nspinhalfs; i++) {
		auto spini = b.spinop(i);
		for(int j = 0; j < i; j++) {
			auto spinj = b.spinop(j);
			res += s.Jin[idx]*(0.5*(spini[0]*spinj[1] + spini[1]*spinj[0]) + spini[2]*spinj[2]);
			idx++;
		}
	}
	return res;
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
	max_worm_count_ = std::max(si.basis.worms.size(), sj.basis.worms.size());

	transitions_.resize(legstates_.size() * max_worm_count_ * leg_count);

	std::vector<int> steps;
	std::vector<int> inv_steps(leg_count*max_worm_count_);
	std::vector<int> step_idx(leg_count*max_worm_count_,-1);

	for(worm_idx worm = 0; worm < max_worm_count_; worm++) {
		for(int leg = 0; leg < leg_count; leg++) {
			const auto &basis = (leg&1 ? sj : si).basis;
			if(worm < static_cast<int>(basis.worms.size())) {
				steps.push_back(worm*leg_count+leg);
				step_idx[worm*leg_count+leg] = steps.size()-1;
				inv_steps[worm*leg_count+leg] = basis.worms[worm].inverse_idx*leg_count+leg;
			}
		}
	}
	struct vertex_change {
		const std::vector<int> &inv_steps;

		int step_in{};
		int step_out{};


		vertex_change inverse() const {
			return vertex_change{inv_steps, inv_steps[step_out], inv_steps[step_in]};
		}
		bool operator==(const vertex_change &other) const {
			return step_in == other.step_in && step_out == other.step_out;
		}
	};

	std::vector<vertex_change> variables;
	for(int step_in : steps) {
		for(int step_out : steps) {
			vertex_change vc{inv_steps, step_in, step_out};

			auto it = std::find_if(variables.begin(), variables.end(),
			                       [&](auto x) { return vc.inverse() == x; });
			if(it == variables.end()) {
				variables.push_back(vc);
			}
		}
	}

	std::vector<solp::constraint> constraints;
	std::vector<double> objective(variables.size(), 0);
	std::transform(
	    variables.begin(), variables.end(), objective.begin(), [&](const vertex_change &vc) {
		bool exclude_bounce = vc.inverse() == vc;
		return exclude_bounce;
	    });

	for(int step_out : steps) {
		std::vector<double> coeff(variables.size(), 0);

		for(size_t i = 0; i < variables.size(); i++) {
			coeff[i] =
			    (variables[i].step_in == step_out) ||
			    (variables[i].inverse().step_in == step_out);
		}
		constraints.push_back(solp::constraint{coeff, 0});
	}

	for(auto &t : transitions_) {
		t.probs.resize(leg_count*max_worm_count_,0);
		t.targets.resize(leg_count*max_worm_count_);
	}

	for(size_t v = 0; v < weights_.size(); v++) {
		for(int step_in : steps) {
			std::vector<int> targets(leg_count*max_worm_count_);

			for(int step_out : steps) {
				int leg_in = step_in % leg_count, leg_out = step_out % leg_count;
				worm_idx worm_in = step_in / leg_count, worm_out = step_out / leg_count;
				int target = vertex_change_apply(si.basis, sj.basis, v, leg_in, worm_in, leg_out, worm_out);

				targets[step_out] = target;
				constraints[step_idx[step_out]].rhs = target >= 0 ? weights_[target] : 0;
			}
			solp::result result = solp::solve(objective, constraints, solp::options{tolerance});

			for(size_t i = 0; i < variables.size(); i++) {
				const auto &var = variables[i];
				int in = var.step_in;
				int out = var.step_out;
				int in_inv = var.inverse().step_in;

				assert(result.x[i] >= -tolerance);
				if(result.x[i] < tolerance) {
					result.x[i] = 0;
				}

				if(targets[in] >= 0) {
					transitions_[targets[in] * max_worm_count_ * leg_count + in].targets[out] =
					    wrap_vertex_idx(targets[in_inv]);
					double norm = constraints[step_idx[in]].rhs == 0 ? 1 : constraints[step_idx[in]].rhs;
					assert(norm > 0);
					transitions_[targets[in] * max_worm_count_ * leg_count + in].probs[out] =
					    result.x[i]/norm;
				}

				in = var.inverse().step_in;
				out = var.inverse().step_out;
				in_inv = var.step_in;

				if(targets[in] >= 0) {
					transitions_[targets[in] * max_worm_count_ * leg_count + in].targets[out] =
					    wrap_vertex_idx(targets[in_inv]);
					double norm = constraints[step_idx[in]].rhs == 0 ? 1 : constraints[step_idx[in]].rhs;
					assert(norm > 0);
					transitions_[targets[in] * max_worm_count_ * leg_count + in].probs[out] =
					    result.x[i]/norm;
				}
			}
		}
	}

	for(auto &t : transitions_) {
		int idx{};
		for(auto &p : t.probs) {
			if(p > 0 && t.targets[idx].invalid()) {
				t.print();
				throw std::runtime_error{fmt::format("it is possible to reach an invalid vertex with p={}", p)};
			}
			idx++;
		}

		std::partial_sum(t.probs.begin(), t.probs.end(), t.probs.begin());
		assert(t.probs.back() < 1.0000000001);
		//assert(!t.invalid());
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
		std::vector<uint32_t> codes;
		for(const auto &t : targets) {
			codes.push_back(t.vertex_idx());
		}
		std::cout << fmt::format("{:.2f}|{:2d}\n",
		                         fmt::join(probs.begin(), probs.end(), " "),
		                         fmt::join(codes.begin(), codes.end(), " "));
}
	

void vertex_data::print(const site_basis &bi, const site_basis &bj) const {
	for(size_t v = 0; v < legstates_.size(); v++) {
		for(int in = 0; in < max_worm_count_ * leg_count; in++) {
			const auto &trans = transitions_[v * max_worm_count_ * leg_count + in];
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
