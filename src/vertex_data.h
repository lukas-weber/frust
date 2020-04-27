#pragma once

#include "bond.h"
#include <vector>
#include <unordered_map>
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
		std::array<int,leg_count*site_basis::worm_count> targets;

		friend vertex_data;
	};

	static constexpr transition invalid_transition{};

	std::vector<std::array<state_idx,4>> legstates;
	double energy_offset{};

	const transition &get_transition(opercode op, int leg_in, worm_idx worm_in) const;
	double get_weight(opercode op) const;
	int get_sign(opercode op) const;
	
	vertex_data(const uc_bond &b, const uc_site &si, const uc_site &sj);
	void print(const site_basis &bi, const site_basis &bj) const;
private:
	
	std::vector<int8_t> signs_;
	std::vector<double> weights_;
	std::vector<transition> transitions_; // [vertex*leg_count*worm_count + worm_in*leg_count + leg_in]
	std::unordered_map<uint32_t, int> code_to_idx_;
	
	void construct_vertices(const uc_bond &b, const uc_site &si, const uc_site &sj, double tolerance);
	void init_code_to_idx();
	int vertex_change_apply(const site_basis &bi, const site_basis &bj, int vertex, int leg_in, worm_idx worm_in_idx, int leg_out, worm_idx worm_out_idx) const;
};

inline int vertex_data::get_sign(opercode op) const {
	int v = code_to_idx_.at(op.vertex());
	return signs_[v];
}

inline const vertex_data::transition &vertex_data::get_transition(opercode op, int leg_in, worm_idx worm_in) const {
	if(code_to_idx_.count(op.vertex()) == 0) {
		return vertex_data::invalid_transition;
	}
		
	int v = code_to_idx_.at(op.vertex());
	return transitions_[v*leg_count*site_basis::worm_count + worm_in*leg_count + leg_in];
}

inline double vertex_data::get_weight(opercode op) const {
	if(code_to_idx_.count(op.vertex()) == 0) {
		return 0;
	}
	int v = code_to_idx_.at(op.vertex());
	return weights_[v];
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

	assert(targets[out] >= 0);

	return std::tuple{leg, worm_out, targets[out]};
}
