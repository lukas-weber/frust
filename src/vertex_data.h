#pragma once

#include "opercode.h"
#include "worms.h"
#include <Eigen/Dense>
#include <vector>

class vertex_data {
public:
	const int leg_count{};

	class transition {
	public:
		bool invalid() const {
			return probs[probs.size() - 1] < 1 - 1e-7;
		}
		void print() const;

	private:
		std::vector<double> probs;
		std::vector<vertexcode> targets;

		friend vertex_data;
	};

	static const transition invalid_transition;

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

	std::vector<state_idx> legstates_; // [vertex; leg_count]

	vertexcode wrap_vertex_idx(int vertex_idx);
	int vertex_change_apply(int vertex, int leg_in, worm_idx worm_in, int leg_out,
	                        worm_idx worm_out) const;
};

inline const vertex_data::transition invalid_transition{};

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
	uint32_t out;
	for(out = 0; out < t.probs.size(); out++) {
		if(random < t.probs[out]) {
			break;
		}
	}

	assert(out < t.probs.size());

	assert(!t.targets[out].invalid());

	if(leg_count == 4) {
		int leg = out % leg_count;
		worm_idx worm_out = out / leg_count;
		return std::tuple{leg, worm_out, t.targets[out]};
	}
	if(leg_count == 6) {
		int leg = out % leg_count;
		worm_idx worm_out = out / leg_count;
		return std::tuple{leg, worm_out, t.targets[out]};
	}
	assert(false);
	__builtin_unreachable();
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
