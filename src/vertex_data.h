#pragma once

#include "bond.h"
#include "opercode.h"
#include "worms.h"
#include <Eigen/Dense>
#include <vector>

class vertex_data {
public:
	static const int leg_count = 4;

	class transition {
	public:
		bool invalid() const {
			return probs[probs.size() - 1] < 1 - 1e-8;
		}
		void print() const;

	private:
		std::vector<double> probs;
		std::vector<vertexcode> targets;

		friend vertex_data;
	};

	static const transition invalid_transition;

	double energy_offset{};

	auto scatter(vertexcode v, int leg_in, worm_idx worm_in, double random) const;
	double get_weight(vertexcode v) const;
	int get_sign(vertexcode v) const;
	vertexcode get_diagonal_vertex(state_idx state_i, state_idx state_j) const;
	const std::array<state_idx, 4> &get_legstate(vertexcode v) const;
	int vertex_count() const;
	
	// bond_hamiltonian has dimension [dim_i;dim_j;dim_i;dim_j]
	vertex_data(int dim_i, int dim_j, const Eigen::MatrixXd& bond_hamiltonian);
	//void print(const site_basis &bi, const site_basis &bj) const;

private:
	int dim_j_{};
	int max_worm_count_{};
	std::vector<vertexcode>
	    diagonal_vertices_; // [dim_i_; dim_j_]
	    
	std::vector<int8_t> signs_;
	std::vector<double> weights_;
	std::vector<transition>
	    transitions_; // [vertex*leg_count*worm_count + worm_in*leg_count + leg_in]

	std::vector<std::array<state_idx, leg_count>> legstates_;

	void construct_vertices(int dim_i, int dim_j, const Eigen::MatrixXd& bond_hamiltonian, double tolerance);
	vertexcode wrap_vertex_idx(int vertex_idx);
	int vertex_change_apply(int dim_i, int dim_j, int vertex, int leg_in,
	                        worm_idx worm_in, int leg_out, worm_idx worm_out) const;
};

inline const vertex_data::transition invalid_transition{};

inline int vertex_data::vertex_count() const {
	return weights_.size();
}

inline vertexcode vertex_data::get_diagonal_vertex(state_idx state_i, state_idx state_j) const {
	return diagonal_vertices_[state_i * dim_j_ + state_j];
}

inline int vertex_data::get_sign(vertexcode v) const {
	return signs_[v.vertex_idx()];
}

inline auto vertex_data::scatter(vertexcode v, int leg_in, worm_idx worm_in, double random) const {
	int vi = v.vertex_idx();
	const auto &t = transitions_[vi * leg_count * max_worm_count_ + worm_in * leg_count + leg_in];

	assert(!t.invalid());
	uint32_t out;
	for(out = 0; out < t.probs.size(); out++) {
		if(random < t.probs[out]) {
			break;
		}
	}

	assert(out < t.probs.size());

	int leg = out % leg_count;
	worm_idx worm_out = out / leg_count;

	assert(!t.targets[out].invalid());

	return std::tuple{leg, worm_out, t.targets[out]};
}

inline double vertex_data::get_weight(vertexcode v) const {
	if(v.invalid()) {
		return 0;
	}
	return weights_[v.vertex_idx()];
}

inline const std::array<state_idx, 4> &vertex_data::get_legstate(vertexcode v) const {
	return legstates_[v.vertex_idx()];
}
