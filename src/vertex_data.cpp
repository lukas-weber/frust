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
using bondmat = Eigen::Matrix<double, jm::basis_size * jm::basis_size, jm::basis_size *jm::basis_size>;

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
	if(nspinhalfs == 1) {
		return mat::Zero();
	} else if(nspinhalfs == 2) {
		auto S = spinops(nspinhalfs);

		mat S2 = 0.5 * (S[0] * S[1] + S[1] * S[0]) + S[2] * S[2];
		return 0.5*(S2 - mat::Identity()*1.5);
	}
	assert(false);
	return mat::Zero();
}

static auto spinhalfops(int site, int nspinhalfs) {
	assert(site >= 0 && site < nspinhalfs);
	auto Shalf = spinops(1);
	if(nspinhalfs == 1) {
		return Shalf;
	}

	double sh = sqrt(0.5);
	// ++ +- -+ --
	Eigen::Matrix<double, jm::basis_size, 4> clebsch_gordan;
	clebsch_gordan << 
		0, sh, -sh, 0,
		1, 0, 0, 0,
		0, sh, sh, 0,
		0, 0, 0, 1;

	decltype(Shalf) S;
	for(size_t i = 0; i < S.size(); i++) {
		S[i] = mat::Zero();
		for(int j1 = 0; j1 < jm::basis_size; j1++) {
			for(int j2 = 0; j2 < jm::basis_size; j2++) {
				for(int m = 0; m < 2; m++) {
					for(int m1 = 0; m1 < 2; m1++) {
						for(int m2 = 0; m2 < 2; m2++) {
							int midx1 = site == 0 ? 2*m1+m : 2*m+m1;
							int midx2 = site == 0 ? 2*m2+m : 2*m+m2;

							S[i](j1,j2) += clebsch_gordan(j1,midx1)*clebsch_gordan(j2,midx2)*Shalf[i](m1,m2);
						}
					}
				}
			}
		}
	}

	return S;
}

static bondmat scalar_product(const std::array<mat,3> &Si, const std::array<mat,3> &Sj) {
	return 0.5 * (kronecker_prod(Si[0], Sj[1]) + kronecker_prod(Si[1], Sj[0])) +
	           kronecker_prod(Si[2], Sj[2]);
}

void vertex_data::construct_vertices(const bond &b, const site &si, const site &sj, double tolerance) {
	std::vector<mat> S(3);

	bondmat H = bondmat::Zero();
	mat Id = mat::Identity();
	H += si.Jin * kronecker_prod(jin_term(si.nspinhalfs),Id);
	H += sj.Jin * kronecker_prod(Id,jin_term(sj.nspinhalfs));

	assert(static_cast<int>(b.J.size()) == si.nspinhalfs*sj.nspinhalfs);
	for(int spini = 0; spini < si.nspinhalfs; spini++) {
		for(int spinj = 0; spinj < sj.nspinhalfs; spinj++) {
			auto Shalfi = spinhalfops(spini, si.nspinhalfs);
			auto Shalfj = spinhalfops(spinj, sj.nspinhalfs);

			H += b.J[sj.nspinhalfs*spini + spinj] * scalar_product(Shalfi, Shalfj);
		}
	}
	
	H.diagonal() *= -1; // from -β ;)

	double epsilon = 0;
	energy_offset = H.diagonal().minCoeff()-epsilon;
	H -= decltype(H)::Identity()*energy_offset;

	for(uint8_t l0 = 0; l0 < jm::basis_size; l0++) {
		for(uint8_t l1 = 0; l1 < jm::basis_size; l1++) {
			for(uint8_t l2 = 0; l2 < jm::basis_size; l2++) {
				for(uint8_t l3 = 0; l3 < jm::basis_size; l3++) {
					double w = H(l0 * jm::basis_size + l1, l2 * jm::basis_size + l3);
					int lbi = jm::local_basis_size(si);
					int lbj = jm::local_basis_size(sj);
					if(fabs(w) > tolerance  && l0 < lbi && l1 < lbj && l2 < lbi && l3 < lbj) {
						legstates.push_back({jm{l0}, jm{l1}, jm{l2}, jm{l3}});
						weights_.push_back(fabs(w));
						signs_.push_back(w >= 0 ? 1 : -1);
					}
				}
			}
		}
	}
}

int vertex_data::vertex_change_apply(const site &si, const site &sj, int vertex, int leg_in, int func_in_idx, int leg_out, int func_out_idx) const {
	auto legstate = legstates.at(vertex);

	const worm_function &func_in = wormfuncs[func_in_idx];
	const worm_function &func_out = wormfuncs[func_out_idx];
	
	const auto &site_in = leg_in&1 ? sj : si;
	const auto &site_out = leg_out&1 ? sj : si;

	legstate[leg_in] = func_in(site_in, legstate[leg_in]);
	if(legstate[leg_in] == jm::invalid) {
		if(leg_in == leg_out && func_in_idx == func_out.inverse_idx) {
			return vertex;
		}
		return -1;
	}
		
	legstate[leg_out] = func_out(site_out, legstate[leg_out]);

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
	const double tolerance = 1e-10;
	construct_vertices(b, si, sj, tolerance);
	transitions_.resize(legstates.size() * wormfunc_count * leg_count);

	struct vertex_change {
		int leg_in{};
		int func_in{};
		int leg_out{};
		int func_out{};

		vertex_change inverse() const {
			return vertex_change{leg_out, wormfuncs[func_out].inverse_idx, leg_in, wormfuncs[func_in].inverse_idx};
		}
		bool operator==(const vertex_change &other) const {
			return leg_in == other.leg_in && func_in == other.func_in &&
			       leg_out == other.leg_out && func_out == other.func_out;
		}
	};

	std::vector<vertex_change> variables;
	for(int func_in = 0; func_in < wormfunc_count; func_in++) {
		for(int leg_in = 0; leg_in < leg_count; leg_in++) {
			for(int func_out = 0; func_out < wormfunc_count; func_out++) {
				for(int leg_out = 0; leg_out < leg_count; leg_out++) {
					vertex_change vc{leg_in, func_in, leg_out, func_out};

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
		    bool exclude_bounce = vc.leg_in == vc.leg_out && vc.func_in == wormfuncs[vc.func_out].inverse_idx;
		    bool exclude_singtrip = wormfuncs[vc.func_out].name[0] == 'D';
		    return exclude_bounce + 0.01*exclude_singtrip;
	    });

	for(int func_out = 0; func_out < wormfunc_count; func_out++) {
		for(int leg_out = 0; leg_out < leg_count; leg_out++) {
			std::vector<double> coeff(variables.size(), 0);

			for(size_t i = 0; i < variables.size(); i++) {
				coeff[i] =
				    (variables[i].leg_in == leg_out && variables[i].func_in == func_out) ||
				    (variables[i].inverse().leg_in == leg_out && variables[i].inverse().func_in == func_out);
			}
			constraints.push_back(solp::constraint{coeff, 0});
		}
	}

	for(auto &t : transitions_) {
		t.probs.fill(-1);
	}

	for(size_t v = 0; v < weights_.size(); v++) {
		for(int func_in = 0; func_in < wormfunc_count; func_in++) {
			for(int leg_in = 0; leg_in < leg_count; leg_in++) {
				std::vector<double> ws(leg_count * wormfunc_count); // [func_out*leg_count+leg_out]
				std::vector<int> targets(leg_count * wormfunc_count);

				for(int func_out = 0; func_out < wormfunc_count; func_out++) {
					for(int leg_out = 0; leg_out < leg_count; leg_out++) {
						int target = vertex_change_apply(si, sj, v, leg_in, func_in, leg_out, func_out);
						int out = func_out * leg_count + leg_out;

						targets[out] = target;
						constraints[out].rhs = target >= 0 ? weights_[target] : 0;
					}
				}
				solp::result result = solp::solve(objective, constraints, solp::options{tolerance});

				for(size_t i = 0; i < variables.size(); i++) {
					const auto &var = variables[i];
					int in = var.func_in * leg_count + var.leg_in;
					int out = var.func_out * leg_count + var.leg_out;
					int in_inv = var.inverse().func_in * leg_count + var.inverse().leg_in;

					assert(result.x[i] >= -tolerance);
					if(result.x[i] < 0) {
						result.x[i] = 0;
					}

					if(targets[in] >= 0) {
						transitions_[targets[in] * wormfunc_count * leg_count + in].targets[out] =
						    targets[in_inv];
						double norm = constraints[in].rhs == 0 ? 1 : constraints[in].rhs;
						assert(norm > 0);
						transitions_[targets[in] * wormfunc_count * leg_count + in].probs[out] =
						    result.x[i]/norm;
					}

					in = var.inverse().func_in * leg_count + var.inverse().leg_in;
					out = var.inverse().func_out * leg_count + var.inverse().leg_out;
					in_inv = var.func_in * leg_count + var.leg_in;

					if(targets[in] >= 0) {
						transitions_[targets[in] * wormfunc_count * leg_count + in].targets[out] =
						    targets[in_inv];
						double norm = constraints[in].rhs == 0 ? 1 : constraints[in].rhs;
						assert(norm > 0);
						transitions_[targets[in] * wormfunc_count * leg_count + in].probs[out] =
						    result.x[i]/norm;
					}
				}
			}
		}
	}

	for(auto &t : transitions_) {
		std::partial_sum(t.probs.begin(), t.probs.end(), t.probs.begin());
		assert(t.probs.back() < 1.0000000001);
		assert(!t.invalid());
	}

	init_code_to_idx();
}

void vertex_data::transition::print() const {
		auto tmp = probs;
		for(size_t i = 1; i < probs.size(); i++) {
			tmp[i] = probs[i]-probs[i-1];
			if(tmp[i] > -1e-8) {
				tmp[i] = fabs(tmp[i]);
			}
		}
		std::cout << fmt::format("{:.2f}|{:2d}\n",
		                         fmt::join(tmp.begin(), tmp.end(), " "),
		                         fmt::join(targets.begin(), targets.end(), " "));
}
	

void vertex_data::print(const site &si, const site &sj) const {
	for(size_t v = 0; v < legstates.size(); v++) {
		for(size_t in = 0; in < wormfunc_count * leg_count; in++) {
			const auto &trans = transitions_[v * wormfunc_count * leg_count + in];
			trans.print();
		}
		std::cout << "\n";
	}
	for(size_t v = 0; v < legstates.size(); v++) {
		for(size_t in = leg_count; in < wormfunc_count * leg_count; in++) {
			int idx = v * wormfunc_count * leg_count +in;
			double p = transitions_[idx].probs[in] - transitions_[idx].probs[in-1];
			if(p > 1e-10) {
				//std::cout << fmt::format("bounce {}-{}: {}\n", v, in, p);
			}
		}
	}

	std::cout << "\nlegstates:\n";
	std::vector<int> idxs(weights_.size());
	int idx{};
	for(const auto &ls : legstates) {
		(void)(ls);
		idxs[idx] = idx;
		idx++;
	}
	std::sort(idxs.begin(), idxs.end(), [&](int i, int j) {return weights_[i] < weights_[j];});


	for(const auto &idx : idxs) {
		const auto &ls = legstates[idx];
		std::cout << fmt::format("{}({}): [{}{} {}{}] ~ {}\n", idx, signs_[idx] > 0 ? '+' : '-', ls[0].name(si), ls[1].name(sj), ls[2].name(si), ls[3].name(sj), weights_[idx]);
	}

}

