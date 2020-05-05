#pragma once

#include "bond.h"
#include <vector>
#include "opercode.h"

class vertex_data {
public: 
	static const int leg_count = 4;
	
	class transition {
	public:
		auto scatter(double randnum) const;
		bool invalid() const {
			return probs[probs.size()-1] < 1-1e-8;
		}
		void print() const;
	private:
		std::array<double,leg_count*site_basis::worm_count> probs;
		std::array<vertexcode,leg_count*site_basis::worm_count> targets;

		friend vertex_data;
	};

	static const transition invalid_transition;

	double energy_offset{};

	const transition &get_transition(vertexcode v, int leg_in, worm_idx worm_in) const;
	double get_weight(vertexcode v) const;
	int get_sign(vertexcode v) const;
	vertexcode get_diagonal_vertex(state_idx state_i, state_idx state_j) const;
	const std::array<state_idx, 4> &get_legstate(vertexcode v) const;
	int vertex_count() const;
	
	vertex_data(const uc_bond &b, const uc_site &si, const uc_site &sj);
	void print(const site_basis &bi, const site_basis &bj) const;
private:
	std::array<vertexcode,site_basis::max_size*site_basis::max_size> diagonal_vertices_; // [site_basis::max_size * statei + state_j]
	
	std::vector<int8_t> signs_;
	std::vector<double> weights_;
	std::vector<transition> transitions_; // [vertex*leg_count*worm_count + worm_in*leg_count + leg_in]
	
	std::vector<std::array<state_idx,4>> legstates_;

	void construct_vertices(const uc_bond &b, const uc_site &si, const uc_site &sj, double tolerance);

	vertexcode wrap_vertex_idx(int vertex_idx);
	int vertex_change_apply(const site_basis &bi, const site_basis &bj, int vertex, int leg_in, worm_idx worm_in_idx, int leg_out, worm_idx worm_out_idx) const;
};

inline const vertex_data::transition invalid_transition{};

inline int vertex_data::vertex_count() const {
	return weights_.size();
}

inline vertexcode vertex_data::get_diagonal_vertex(state_idx state_i, state_idx state_j) const {
	return diagonal_vertices_[state_i*site_basis::max_size + state_j];
}

inline int vertex_data::get_sign(vertexcode v) const {
	return signs_[v.vertex_idx()];
}

inline const vertex_data::transition &vertex_data::get_transition(vertexcode v, int leg_in, worm_idx worm_in) const {
	int vi = v.vertex_idx();
	return transitions_[vi*leg_count*site_basis::worm_count + worm_in*leg_count + leg_in];
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

inline auto vertex_data::transition::scatter(double randnum) const {
	assert(!invalid());
	uint32_t out;
	for(out = 0; out < probs.size(); out++) {
		if(randnum < probs[out]) {
			break;
		}
	}

	assert(out < leg_count * site_basis::worm_count);

	int leg = out%leg_count;
	worm_idx worm_out = out/leg_count;

	assert(!targets[out].invalid());

	return std::tuple{leg, worm_out, targets[out]};
}
