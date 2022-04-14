#include "vertex_data.h"

#include <Highs.h>
#include <fmt/format.h>
#include <iostream>
#include <numeric>

static double calc_energy_offset(const Eigen::MatrixXd &H) {
	double hmin = H.diagonal().minCoeff();
	double hmax = std::max(fabs(hmin), fabs(H.maxCoeff()));

	double epsilon = hmax / 2;
	return hmin - epsilon;
}

static auto construct_vertices(const std::vector<int> &dims, const Eigen::MatrixXd &H,
                               double energy_offset, double tolerance) {
	std::vector<vertexcode> diagonal_vertices(H.cols());
	std::vector<double> weights;
	std::vector<state_idx> legstates;
	std::vector<int8_t> signs;

	for(int i = 0; i < H.rows(); i++) {
		for(int j = 0; j < H.cols(); j++) {
			double w = -H(i, j);
			if(i == j) {
				w -= energy_offset;
			}

			if(fabs(w) > tolerance) {
				if(i == j) {
					diagonal_vertices[i] = vertexcode{true, static_cast<uint32_t>(weights.size())};
				}

				for(int tmp : {i, j}) {
					std::vector<state_idx> legstate;

					for(int k = 0; k < static_cast<int>(dims.size()); k++) {
						int d = dims[dims.size() - k - 1];
						legstate.push_back(tmp % d);
						tmp /= d;
					}
					std::reverse(legstate.begin(), legstate.end());
					std::copy(legstate.begin(), legstate.end(), std::back_inserter(legstates));
				}

				weights.push_back(fabs(w));
				signs.push_back(w >= 0 ? 1 : -1);
			}
		}
	}

	return std::make_tuple(diagonal_vertices, weights, legstates, signs);
}

vertexcode vertex_data::wrap_vertex_idx(int vertex_idx) {
	if(vertex_idx < 0) {
		return vertexcode{};
	}

	const state_idx *ls = &legstates_[leg_count * vertex_idx];
	bool diagonal = true;
	for(int i = 0; i < leg_count / 2; i++) {
		diagonal &= ls[i] == ls[leg_count / 2 + i];
	}
	return vertexcode{diagonal, static_cast<uint32_t>(vertex_idx)};
}

int vertex_data::vertex_change_apply(int vertex, int leg_in, worm_idx worm_in, int leg_out,
                                     worm_idx worm_out) const {
	std::vector<state_idx> new_legstate(&legstates_[leg_count * vertex],
	                                    &legstates_[leg_count * vertex] + leg_count);

	int dim_in = dims[leg_in % (leg_count / 2)];
	int dim_out = dims[leg_out % (leg_count / 2)];

	new_legstate[leg_in] = worm_action(worm_in, new_legstate[leg_in], dim_in);
	new_legstate[leg_out] = worm_action(worm_out, new_legstate[leg_out], dim_out);
	for(int v = 0; v < static_cast<int>(weights_.size()); v++) {
		bool match = true;
		for(int l = 0; l < leg_count; l++) {
			match &= new_legstate[l] == legstates_[leg_count * v + l];
		}
		if(match) {
			return v;
		}
	}
	return -1;
}

vertex_data::vertex_data(const std::vector<int> &dims, const Eigen::MatrixXd &bond_hamiltonian)
    : leg_count{2 * static_cast<int>(dims.size())}, dims{dims},
      max_worm_count_{worm_count(*std::max_element(dims.begin(), dims.end()))} {
	assert(leg_count >= 2);
	int total_dim =
	    std::accumulate(dims.begin(), dims.end(), 1, [](int a, int b) { return a * b; });
	assert(total_dim == bond_hamiltonian.cols() && total_dim == bond_hamiltonian.rows());

	const double tolerance = 1e-7;
	energy_offset = calc_energy_offset(bond_hamiltonian);
	std::tie(diagonal_vertices_, weights_, legstates_, signs_) =
	    construct_vertices(dims, bond_hamiltonian, energy_offset, tolerance);

	assert(legstates_.size() == leg_count * weights_.size());
	transitions_.resize(weights_.size() * max_worm_count_ * leg_count);

	std::vector<int> steps;
	std::vector<int> inv_steps(leg_count * max_worm_count_);
	std::vector<int> step_idx(leg_count * max_worm_count_, -1);

	for(worm_idx worm = 0; worm < max_worm_count_; worm++) {
		for(int leg = 0; leg < leg_count; leg++) {
			int dim = dims[leg % (leg_count / 2)];
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
	model.lp_.a_matrix_.value_ = std::vector<double>(constraints_col_idxs.size(), 1);

	for(auto &t : transitions_) {
		t.probs.resize(leg_count * max_worm_count_, 0);
		t.targets.resize(leg_count * max_worm_count_);
	}

	std::vector<double> constraints(steps.size());

	Highs highs;
	highs.setOptionValue("log_dev_level", kHighsLogDevLevelNone);
	highs.setOptionValue("primal_feasibility_tolerance", 1e-10);
	highs.setOptionValue("dual_feasibility_tolerance", 1e-10);
	highs.setOptionValue("output_flag", false);
	for(size_t v = 0; v < weights_.size(); v++) {
		for(int step_in : steps) {
			std::vector<int> targets(leg_count * max_worm_count_);

			for(int step_out : steps) {
				int leg_in = step_in % leg_count, leg_out = step_out % leg_count;
				worm_idx worm_in = step_in / leg_count, worm_out = step_out / leg_count;
				int target = vertex_change_apply(v, leg_in, worm_in, leg_out, worm_out);

				targets[step_out] = target;
				constraints[step_idx[step_out]] = target >= 0 ? weights_[target] : 0;
			}

			model.lp_.row_lower_ = constraints;
			model.lp_.row_upper_ = constraints;

			HighsStatus status = highs.passModel(model);
			assert(status == HighsStatus::kOk);
			status = highs.run();
			assert(status == HighsStatus::kOk);

			const auto &model_status = highs.getModelStatus();
			assert(model_status == HighsModelStatus::kOptimal);
			const auto &solution = highs.getSolution();

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
					double norm = constraints[step_idx[in]] == 0 ? 1 : constraints[step_idx[in]];
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
					double norm = constraints[step_idx[in]] == 0 ? 1 : constraints[step_idx[in]];
					assert(norm > 0);
					transitions_[targets[in] * max_worm_count_ * leg_count + in].probs[out] =
					    prob / norm;
				}
			}
		}
	}

	for(auto &t : transitions_) {
		int idx{};
		double norm{};
		for(auto &p : t.probs) {
			if(p > 0 && t.targets[idx].invalid()) {
				t.print();
				throw std::runtime_error{
				    fmt::format("it is possible to reach an invalid vertex with p={}", p)};
			}
			norm += p;
			idx++;
		}
		if(norm > 0) {
			assert(fabs(norm - 1) <
			       1e-7); // as long as the norm is not completely off, renormalize.
			for(auto &p : t.probs) {
				p /= norm;
			}
		}

		std::partial_sum(t.probs.begin(), t.probs.end(), t.probs.begin());
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

void vertex_data::print() const {
	for(size_t v = 0; v < weights_.size(); v++) {
		for(int in = 0; in < max_worm_count_ * leg_count; in++) {
			const auto &trans = transitions_[v * max_worm_count_ * leg_count + in];
			trans.print();
		}
		std::cout << "\n";
	}

	std::cout << "\nlegstates:\n";
	std::vector<int> idxs(weights_.size());
	int idx{};
	for(const auto &w : weights_) {
		(void)(w);
		idxs[idx] = idx;
		idx++;
	}
	// std::sort(idxs.begin(), idxs.end(), [&](int i, int j) {return weights_[i] < weights_[j];});

	for(const auto &idx : idxs) {
		const auto &ls = &legstates_[leg_count * idx];
		std::cout << fmt::format("{}({}): [{}{} {}{}] ~ {}\n", idx, signs_[idx] > 0 ? '+' : '-',
		                         ls[0], ls[1], ls[2], ls[3], weights_[idx]);
	}
}
