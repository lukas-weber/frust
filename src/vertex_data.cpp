#include "vertex_data.h"

#include <fmt/format.h>
#include <iostream>
#include <numeric>
#include <Highs.h>

static double calc_energy_offset(const Eigen::MatrixXd &H) {
	double hmax = 0;
	double hmin = 0;

	
	for(double h : H.diagonal()) {
		if(-h < hmin) {
			hmin = -h;
		}
	}

	for(double h : H.reshaped()) {
		if(fabs(h) > hmax) {
			hmax = h;
		}
	}

	double epsilon = hmax / 2;
	return hmin - epsilon;
}

void vertex_data::construct_vertices(int dim_i, int dim_j, const Eigen::MatrixXd& H,
                                     double tolerance) {
	for(int s = 0; s < H.cols()*H.rows(); s++) {
		state_idx l0 = s / (dim_j * dim_i * dim_j);
		state_idx l1 = (s / (dim_i * dim_j)) % dim_j; 
		state_idx l2 = (s / dim_j) % dim_i;
		state_idx l3 = s % dim_j;

		double w = -H(l0*dim_j + l1, l2*dim_j  + l3);
		if(l0 == l2 && l1 == l3) {
			w -= energy_offset;
		}

		if(fabs(w) > tolerance) {
			if(l0 == l2 && l1 == l3) {
				diagonal_vertices_[l0 * dim_j + l1] = vertexcode{true, static_cast<uint32_t>(weights_.size())};
			}
			legstates_.push_back({l0, l1, l2, l3});
			weights_.push_back(fabs(w));
			signs_.push_back(w >= 0 ? 1 : -1);
		}		
	}
}

vertexcode vertex_data::wrap_vertex_idx(int vertex_idx) {
	if(vertex_idx < 0) {
		return vertexcode{};
	}

	const auto &ls = legstates_[vertex_idx];
	return vertexcode{ls[0] == ls[2] && ls[1] == ls[3], static_cast<uint32_t>(vertex_idx)};
}

int vertex_data::vertex_change_apply(int dim_i, int dim_j, int vertex,
                                     int leg_in, worm_idx worm_in, int leg_out,
                                     worm_idx worm_out) const {
	auto legstate = legstates_.at(vertex);

	int dim_in = leg_in & 1 ? dim_j : dim_i;
	int dim_out = leg_out & 1 ? dim_j : dim_i;

	legstate[leg_in] = worm_action(worm_in, legstate[leg_in], dim_in);
	legstate[leg_out] = worm_action(worm_out, legstate[leg_out], dim_out);
	auto it = std::find(legstates_.begin(), legstates_.end(), legstate);
	if(it == legstates_.end()) {
		return -1;
	}

	return it - legstates_.begin();
}

vertex_data::vertex_data(int dim_i, int dim_j, const Eigen::MatrixXd& bond_hamiltonian) 
	: dim_j_{dim_j}, max_worm_count_{std::max(worm_count(dim_i), worm_count(dim_j))}, diagonal_vertices_(dim_i*dim_j) {
	const double tolerance = 1e-10;
	energy_offset = calc_energy_offset(bond_hamiltonian);
	construct_vertices(dim_i, dim_j, bond_hamiltonian, tolerance);

	transitions_.resize(legstates_.size() * max_worm_count_ * leg_count);

	std::vector<int> steps;
	std::vector<int> inv_steps(leg_count * max_worm_count_);
	std::vector<int> step_idx(leg_count * max_worm_count_, -1);

	for(worm_idx worm = 0; worm < max_worm_count_; worm++) {
		for(int leg = 0; leg < leg_count; leg++) {
			int dim = leg & 1 ? dim_j : dim_i;
			if(worm < worm_count(dim)) {
				steps.push_back(worm * leg_count + leg);
				step_idx[worm * leg_count + leg] = steps.size() - 1;
				inv_steps[worm * leg_count + leg] = worm_inverse(worm, dim) * leg_count + leg;
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

	HighsModel model;
	model.lp_.num_col_ = variables.size();
	model.lp_.num_row_ = steps.size();
	model.lp_.sense_ = ObjSense::kMinimize;
	model.lp_.offset_ = 0;

	std::vector<double> objective(variables.size(), 0);
	std::transform(variables.begin(), variables.end(), objective.begin(),
	               [&](const vertex_change &vc) {
		               bool exclude_bounce = vc.inverse() == vc;
		               return exclude_bounce;
	               });
	model.lp_.col_cost_ = objective;
	model.lp_.col_lower_ = std::vector<double>(variables.size(), 0);
	model.lp_.col_upper_ = std::vector<double>(variables.size(), 1e30);
	std::vector<int> constraints_row_starts = {0};
	std::vector<int> constraints_col_idxs;
	for(int step_out : steps) {
		for(size_t i = 0; i < variables.size(); i++) {
			if(variables[i].step_in == step_out || variables[i].inverse().step_in == step_out) {
				constraints_col_idxs.push_back(i);
			}
		}
		constraints_row_starts.push_back(constraints_col_idxs.size());
	}
	model.lp_.a_matrix_.format_ = MatrixFormat::kRowwise;
	model.lp_.a_matrix_.start_ = constraints_row_starts;
	model.lp_.a_matrix_.index_ = constraints_col_idxs;
	model.lp_.a_matrix_.value_ = std::vector<double>(constraints_col_idxs.size(),1);


	for(auto &t : transitions_) {
		t.probs.resize(leg_count * max_worm_count_, 0);
		t.targets.resize(leg_count * max_worm_count_);
	}
	
	std::vector<double> constraints_upper(steps.size());

	Highs highs;
	highs.setOptionValue("log_dev_level", kHighsLogDevLevelNone);
	highs.setOptionValue("output_flag", false);
	for(size_t v = 0; v < weights_.size(); v++) {
		for(int step_in : steps) {
			std::vector<int> targets(leg_count * max_worm_count_);

			for(int step_out : steps) {
				int leg_in = step_in % leg_count, leg_out = step_out % leg_count;
				worm_idx worm_in = step_in / leg_count, worm_out = step_out / leg_count;
				int target =
				    vertex_change_apply(dim_i, dim_j, v, leg_in, worm_in, leg_out, worm_out);

				targets[step_out] = target;
				constraints_upper[step_idx[step_out]] = target >= 0 ? weights_[target] : 0;
			}

			model.lp_.row_lower_ = constraints_upper;
			model.lp_.row_upper_ = constraints_upper;

			HighsStatus status = highs.passModel(model);
			assert(status == HighsStatus::kOk);
			status = highs.run();
			assert(status == HighsStatus::kOk);

			const auto& model_status = highs.getModelStatus();
			assert(model_status == HighsModelStatus::kOptimal);
			const auto& solution = highs.getSolution();
			
			for(size_t i = 0; i < variables.size(); i++) {
				const auto &var = variables[i];
				int in = var.step_in;
				int out = var.step_out;
				int in_inv = var.inverse().step_in;

				double prob = solution.col_value[i];

				assert(prob >= -tolerance);
				if(prob < tolerance) {
					prob = 0;
				}

				if(targets[in] >= 0) {
					transitions_[targets[in] * max_worm_count_ * leg_count + in].targets[out] =
					    wrap_vertex_idx(targets[in_inv]);
					double norm =
					    constraints_upper[step_idx[in]] == 0 ? 1 : constraints_upper[step_idx[in]];
					assert(norm > 0);
					transitions_[targets[in] * max_worm_count_ * leg_count + in].probs[out] =
					    prob / norm;
				}

				in = var.inverse().step_in;
				out = var.inverse().step_out;
				in_inv = var.step_in;

				if(targets[in] >= 0) {
					transitions_[targets[in] * max_worm_count_ * leg_count + in].targets[out] =
					    wrap_vertex_idx(targets[in_inv]);
					double norm =
					    constraints_upper[step_idx[in]] == 0 ? 1 : constraints_upper[step_idx[in]];
					assert(norm > 0);
					transitions_[targets[in] * max_worm_count_ * leg_count + in].probs[out] =
					    prob / norm;
				}
			}
		}
	}

	for(auto &t : transitions_) {
		int idx{};
		for(auto &p : t.probs) {
			if(p > 0 && t.targets[idx].invalid()) {
				t.print();
				throw std::runtime_error{
				    fmt::format("it is possible to reach an invalid vertex with p={}", p)};
			}
			idx++;
		}

		std::partial_sum(t.probs.begin(), t.probs.end(), t.probs.begin());
		assert(t.probs.back() < 1.0000000001);
		// assert(!t.invalid());
	}
}

void vertex_data::transition::print() const {
	auto tmp = probs;
	for(size_t i = 1; i < probs.size(); i++) {
		tmp[i] = probs[i] - probs[i - 1];
		if(tmp[i] > -1e-8) {
			tmp[i] = fabs(tmp[i]);
		}
	}
	std::vector<int> codes;
	for(const auto &t : targets) {
		codes.push_back(t.invalid() ? -1 : t.vertex_idx());
	}
	std::cout << fmt::format("{:.2f}|{:2d}\n", fmt::join(probs.begin(), probs.end(), " "),
	                         fmt::join(codes.begin(), codes.end(), " "));
}
/*
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
	// std::sort(idxs.begin(), idxs.end(), [&](int i, int j) {return weights_[i] < weights_[j];});

	for(const auto &idx : idxs) {
		const auto &ls = legstates_[idx];
		std::cout << fmt::format("{}({}): [{}{} {}{}] ~ {}\n", idx, signs_[idx] > 0 ? '+' : '-',
		                         bi.states[ls[0]].name, bj.states[ls[1]].name,
		                         bi.states[ls[2]].name, bj.states[ls[3]].name, weights_[idx]);
	}
}*/
