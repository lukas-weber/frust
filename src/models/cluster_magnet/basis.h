#pragma once
#include "util/kronecker_product.h"
#include <Eigen/Dense>
#include <array>
#include <cstdint>
#include <vector>

#include "worms.h"

class site_basis {
public:
	struct state {
		char name{};
		double j{};
		double m{};
		double jdim{};
	};

	static const uint32_t state_bits = 3; // max number of bits the state occupies
	static const uint32_t max_size = 8;
	static const state_idx invalid = max_size;

	int nspinhalfs{};
	std::vector<state> states;

	site_basis(int nspinhalfs, const std::vector<state> &states, const Eigen::MatrixXd &trans);

	int size() const;
	double m(state_idx st) const;
	double j(state_idx st) const;
	char name(state_idx st) const;

	std::array<Eigen::MatrixXd, 2> spinop(int spinhalf) const;

private:
	Eigen::MatrixXd trans_; // trans.col(st) is the sz-basis vector of state st.
};

namespace site_bases {

// clang-format off
	const state_idx in_ = site_basis::invalid;
	const double isq_ = 1 / sqrt(2.);
	const double isq3_ = 1 / sqrt(3.);
	const double isq6_ = 1 / sqrt(6.);
	const site_basis spin = {
		1,
		{{'+', 0.5, 0.5}, {'-', 0.5, -0.5}},
		init_mat(2,2, {1,0,0,1}),
	};
	const site_basis dimer = {
		2,
		{{'*', 0, 0}, {'P', 1, 1}, {'0', 1, 0}, {'M', 1, -1}},
		init_mat(4,4, {0,     1, 0,    0,
			       isq_,  0, isq_, 0,
			       -isq_, 0, isq_, 0,
			        0,    0, 0,    1}),
	};
	const site_basis trimer = {
		3,
		{{'a', 0.5, 0.5, 0}, {'b', 0.5, -0.5, 0},
		 {'x', 0.5, 0.5, 1}, {'y', 0.5, -0.5, 1},
		 {'P', 1.5, 1.5, 1}, {'p', 1.5, 0.5, 1}, {'m', 1.5, -0.5, 1}, {'M', 1.5, -1.5, 1}},
		init_mat(8,8, {0,	 0,  isq_,     0, -isq_,     0,	       0, 0,
			       0,	 0,     0,  -isq_,     0, isq_,	       0, 0,
			       0, -2*isq6_, isq6_,     0, isq6_,     0,	       0, 0,
			       0,	 0,	0, isq6_,     0, isq6_, -2*isq6_, 0,
			       1,	 0,	0,     0,     0,     0,	       0, 0,
			       0,    isq3_, isq3_,     0, isq3_,     0,	       0, 0,
			       0,	 0,	0, isq3_,     0, isq3_,	   isq3_, 0,
			       0,	 0,     0,     0,     0,     0,	       0, 1}).transpose(),
	};
	const site_basis trimer23 = {
		3,
		{{'a', 0.5, 0.5, 0}, {'b', 0.5, -0.5, 0},
		 {'x', 0.5, 0.5, 1}, {'y', 0.5, -0.5, 1},
		 {'P', 1.5, 1.5, 1}, {'p', 1.5, 0.5, 1}, {'m', 1.5, -0.5, 1}, {'M', 1.5, -1.5, 1}},
		init_mat(8,8, {0, -isq_,  isq_,  0,	   0,        0,     0,     0,
			       0, 0,      0,     0,        0,  	     isq_,  -isq_, 0,
			       0, isq6_,  isq6_, 0,	   -2*isq6_, 0,     0,     0,
			       0, 0,	  0, 	 -2*isq6_, 0,  	     isq6_, isq6_, 0,
			       1, 0,	  0,     0,        0,        0,     0,     0,
			       0, isq3_,  isq3_, 0,	   isq3_,    0,     0,     0,
			       0, 0,	  0, 	 isq3_,    0,  	     isq3_, isq3_, 0,
			       0, 0,      0,     0,        0,        0,     0,     1}).transpose(),
	};
	
	/*const site_basis spin1 = {
		1,
		{{'+', 1, -1}, {'0', 1, 0}, {'-', 1, -1}},
		init_mat(3,3, {1,0,0,
			       0,1,0,
			       0,0,1}),
		generate_addmod_worms(3),
	};*/
// clang-format on
}

inline int site_basis::size() const {
	return states.size();
}
