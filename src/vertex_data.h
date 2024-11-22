#pragma once

#include "bond.h"
#include <vector>
#include <unordered_map>
#include "jm.h"
#include "opercode.h"

class vertex_data {
public: 
	static const int leg_count = 4;
	
	class transition {
	public:
		auto scatter(double randnum) const;
		bool invalid() const {
			return probs[probs.size()-1] < 0.9;
		}
	private:
		std::array<double,leg_count*jm::basis_size> probs;
		std::array<int,leg_count*jm::basis_size> targets;

		friend vertex_data;
	};

	static constexpr transition invalid_transition{};

	std::vector<std::array<jm,4>> legstates;

	const transition &get_transition(opercode op, int leg_in, jm_action action_in) const;
	double get_weight(opercode op) const;
	
	vertex_data(const bond &b, const site &si, const site &sj);
	void print(const site &si, const site &sj) const;
private:
	double energy_offset_{};

	std::vector<double> weights_;
	std::vector<transition> transitions_; // [vertex*leg_count*basis_size + worm_action_in*leg_count + leg_in]
	std::unordered_map<uint32_t, int> code_to_idx_;
	
	void construct_vertices(const bond &b, const site &si, const site &sj);
	void init_code_to_idx();
	int vertex_change_apply(int vertex, int leg_in, jm_action action_in, int leg_out, jm_action action_out) const;
};

inline const vertex_data::transition &vertex_data::get_transition(opercode op, int leg_in, jm_action action_in) const {
	if(code_to_idx_.count(op.vertex()) == 0) {
		return vertex_data::invalid_transition;
	}
		
	int v = code_to_idx_.at(op.vertex());
	return transitions_[v*leg_count*jm::basis_size + action_in*leg_count + leg_in];
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

	assert(out < leg_count * jm::basis_size);

	int leg = out%leg_count;
	jm_action action_out = out/leg_count;

	return std::tuple{leg, action_out, targets[out]};
}
