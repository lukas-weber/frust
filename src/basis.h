#pragma once
#include <cstdint>
#include <vector>
#include <array>
#include <Eigen/Dense>
#include "util.h"

using state_idx = uint8_t;
using worm_idx = int;

struct worm_function {
	const char *name;
	int inverse_idx;
	std::vector<state_idx> action;

	state_idx operator()(state_idx state) const {
		return action[state];
	}
};

class site_basis {
public:
	struct state {
		char name{};
		double j{};
		double m{};
	};

	static const worm_idx worm_count = 3;
	static const uint32_t state_bits = 3; // max number of bits the state occupies
	static const state_idx invalid = 1<<state_bits;

	int nspinhalfs{};
	std::vector<state> states;
	std::array<worm_function, worm_count> worms;

	site_basis(int nspinhalfs, const std::vector<state> &states, const std::array<worm_function, worm_count> &worms, const Eigen::MatrixXd &trans) : nspinhalfs{nspinhalfs}, states{states}, worms{worms}, trans_{trans} {
	}

	int size() const;
	double m(state_idx st) const;
	double j(state_idx st) const;
	char name(state_idx st) const;

	std::array<Eigen::MatrixXd,3> spinop(int spinhalf) const;
private:
	Eigen::MatrixXd trans_; // trans.col(st) is the sz-basis vector of state st.
};

namespace site_bases {
	const state_idx in_ = site_basis::invalid;
	const double isq_ = 1/sqrt(2.);
	
	const site_basis spin = {
		1,
		{{'+', 0.5, 0.5}, {'-', 0.5, -0.5}},
		{worm_function{"S+", 1, {in_, 0}},
		 worm_function{"S-", 0, {1, in_}},
		 worm_function{"NA", 2, {in_, in_}}},
		init_mat(2,2, {1,0,0,1})
	};
	const site_basis dimer = {
		2,
		{{'*', 0, 0}, {'P', 1, 1}, {'0', 1, 0}, {'M', 1, -1}},
		{worm_function{"S+", 1, {in_, in_, 1, 2}},
		 worm_function{"S-", 0, {in_, 2, 3, in_}},
		 worm_function{"Dz", 2, {2, in_, 0, in_}}},
		init_mat(4,4, {0,     1, 0,    0,
			       isq_,  0, isq_, 0,
			       -isq_, 0, isq_, 0,
			        0,    0, 0,    1})
	};
	/*const site_basis trimer = {
		3,
		{{'a', 0.5, 0.5}, {'b', 0.5, -0.5},
		 {'x', 0.5, 0.5}, {'y', 0.5, -0.5},
		 {'P', 1.5, 1.5}, {'p', 1.5, 0.5}, {'m', 1.5, -0.5}, {'M', 1.5, -1.5}},
		{worm_function{"S+", 1, {in_, in_, 1, 2}},
		 worm_function{"S-", 0, {in_, 2, 3, in_}},
		 worm_function{"Dz", 2, {2, in_, 0, in_}}},
		init_mat(2,2, {1,0,0,1})
	};*/
}
	


inline int site_basis::size() const {
	return states.size();
}
