#pragma once

#include "bond.h"
#include <vector>
#include <unordered_map>
#include "jm.h"
#include "opercode.h"

struct worm_function {
	const char *name;
	int inverse_idx;
	std::array<jm, 2> action_shalf;
	std::array<jm, 4> action;

	jm operator()(const site &site, jm state) const {
		if(site.nspinhalfs == 1) {
			return jm{action_shalf[state.code()]};
		}
		assert(site.nspinhalfs == 2);
		return jm{action[state.code()]};
	}
};

class vertex_data {
public: 
	static const int leg_count = 4;
	static const int wormfunc_count = 3;
	static const std::array<worm_function, wormfunc_count> wormfuncs;
	
	class transition {
	public:
		auto scatter(double randnum) const;
		bool invalid() const {
			return probs[probs.size()-1] < 1-1e-8;
		}
		void print() const;
	private:
		std::array<double,leg_count*wormfunc_count> probs;
		std::array<int,leg_count*wormfunc_count> targets;

		friend vertex_data;
	};

	static constexpr transition invalid_transition{};

	std::vector<std::array<jm,4>> legstates;
	double energy_offset{};

	const transition &get_transition(opercode op, int leg_in, jm_action action_in) const;
	double get_weight(opercode op) const;
	int get_sign(opercode op) const;
	
	vertex_data(const bond &b, const site &si, const site &sj);
	void print(const site &si, const site &sj) const;
private:
	
	std::vector<int8_t> signs_;
	std::vector<double> weights_;
	std::vector<transition> transitions_; // [vertex*leg_count*wormfunc_count + wormfunc_in*leg_count + leg_in]
	std::unordered_map<uint32_t, int> code_to_idx_;
	
	void construct_vertices(const bond &b, const site &si, const site &sj, double tolerance);
	void init_code_to_idx();
	int vertex_change_apply(const site &si, const site &sj, int vertex, int leg_in, int func_in_idx, int leg_out, int func_out_idx) const;
};


// these are defined in some ugly enlarged basis to make them e.g. invertible
inline const std::array<worm_function, vertex_data::wormfunc_count> vertex_data::wormfuncs = {{
	{"+", 1, {jm::invalid,0}, {jm::invalid, jm::invalid, 1, 2}},
	{"-", 0, {1,jm::invalid}, {jm::invalid, 2, 3, jm::invalid}},
	{"D0", 2, {0, 1}, {2, jm::invalid, 0, jm::invalid}}
//	{"D+", 3, {1,jm::invalid}, {1, jm::invalid, jm::invalid, 0}},
//	{"D-", 2, {jm::invalid, 0}, {3, 0, jm::invalid, jm::invalid}}
}};
/*	{"^1", 0, {1, 0, 4, 4}, {1,0,3,2,4,4,4,4,4,4}},
	{"^2", 1, {2, 3, 0, 1}, {2,3,0,1,4,4,4,4,4,4}},
}};*/


inline int vertex_data::get_sign(opercode op) const {
	int v = code_to_idx_.at(op.vertex());
	return signs_[v];
}

inline const vertex_data::transition &vertex_data::get_transition(opercode op, int leg_in, jm_action action_in) const {
	if(code_to_idx_.count(op.vertex()) == 0) {
		return vertex_data::invalid_transition;
	}
		
	int v = code_to_idx_.at(op.vertex());
	return transitions_[v*leg_count*wormfunc_count + action_in*leg_count + leg_in];
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

	assert(out < leg_count * vertex_data::wormfunc_count);

	int leg = out%leg_count;
	int wormfunc_out = out/leg_count;

	assert(targets[out] >= 0);

	return std::tuple{leg, wormfunc_out, targets[out]};
}
