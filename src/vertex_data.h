#pragma once

#include "opercode.h"
#include "worms.h"
#include <Eigen/Dense>
#include <vector>

class vertex_data {
public:
	const int leg_count{};

	struct transition {
		bool invalid() const {
			return offset < 0;
		};

		int offset{-1};
		int length{0};
	};

	double energy_offset{};
	const std::vector<int> dims;

	auto scatter(vertexcode v, int leg_in, worm_idx worm_in, double random) const;

	double get_weight(vertexcode v) const;
	int get_sign(vertexcode v) const;
	vertexcode get_diagonal_vertex(
	    int compound_state_idx) const; // compound_state_idx fits right into diagonal_vertices_
	const state_idx *
	    get_legstate(vertexcode v) const; // [leg_count] XXX: replace by std::span or something
	int vertex_count() const;

	// bond_hamiltonian has dimension [dim_i, dim_j, ...; dim_i, dim_j, ...]
	vertex_data(const std::vector<int> &dims, const Eigen::MatrixXd &bond_hamiltonian);
	void print() const;

private:
	int max_worm_count_{};
	std::vector<vertexcode> diagonal_vertices_; // [dims_[0]; dims_[1]; ...]

	std::vector<int8_t> signs_;
	std::vector<double> weights_;

	std::vector<transition>
	    transitions_; // [vertex*leg_count*worm_count + worm_in*leg_count + leg_in]
	// the transition probabilities are stored in a dense array indexed by transition::offset
	std::vector<double> transition_cumprobs_;
	std::vector<vertexcode> transition_targets_;
	std::vector<std::pair<int, int>> transition_step_outs_; // (leg, worm)

	std::vector<state_idx> legstates_; // [vertex; leg_count]

	vertexcode wrap_vertex_idx(int vertex_idx);
	int vertex_change_apply(int vertex, int leg_in, worm_idx worm_in, int leg_out,
	                        worm_idx worm_out) const;
};

inline int vertex_data::vertex_count() const {
	return weights_.size();
}

inline vertexcode vertex_data::get_diagonal_vertex(int compound_state_idx) const {
	return diagonal_vertices_[compound_state_idx];
}

inline int vertex_data::get_sign(vertexcode v) const {
	return signs_[v.vertex_idx()];
}

inline auto vertex_data::scatter(vertexcode v, int leg_in, worm_idx worm_in, double random) const {
	int vi = v.vertex_idx();
	const auto &t = transitions_[vi * leg_count * max_worm_count_ + worm_in * leg_count + leg_in];

	assert(worm_in < worm_count(dims[leg_in % (leg_count / 2)]));

	assert(!t.invalid());
	int out;
	for(out = 0; out < t.length; out++) {
		if(random < transition_cumprobs_[t.offset + out]) {
			break;
		}
	}

	assert(out < t.length);

	auto [leg_out, worm_out] = transition_step_outs_[t.offset + out];
	return std::tuple{leg_out, worm_out, transition_targets_[t.offset + out]};
}

inline double vertex_data::get_weight(vertexcode v) const {
	if(v.invalid()) {
		return 0;
	}
	return weights_[v.vertex_idx()];
}

inline const state_idx *vertex_data::get_legstate(vertexcode v) const {
	return &legstates_[leg_count * v.vertex_idx()];
}
