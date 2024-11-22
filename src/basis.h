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
};

std::vector<worm_function> generate_xor_worms(uint32_t basis_size);
std::vector<worm_function> generate_addmod_worms(uint32_t basis_size);

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
	std::vector<worm_function> worms;

	site_basis(int nspinhalfs, const std::vector<state> &states, const Eigen::MatrixXd &trans , const std::vector<worm_function> &worms);

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
	const double isq3_ = 1/sqrt(3.);
	const double isq6_ = 1/sqrt(6.);

	const site_basis spin = {
		1,
		{{'+', 0.5, 0.5}, {'-', 0.5, -0.5}},
		init_mat(2,2, {1,0,0,1}),
		generate_xor_worms(1<<1),
	};
	const site_basis dimer = {
		2,
		{{'*', 0, 0}, {'P', 1, 1}, {'0', 1, 0}, {'M', 1, -1}},
		init_mat(4,4, {0,     1, 0,    0,
			       isq_,  0, isq_, 0,
			       -isq_, 0, isq_, 0,
			        0,    0, 0,    1}),
		generate_xor_worms(1<<2)
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
		generate_xor_worms(1<<3)       
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
		generate_xor_worms(1<<3)       
	};
	
	/*const site_basis spin1 = {
		1,
		{{'+', 1, -1}, {'0', 1, 0}, {'-', 1, -1}},
		init_mat(3,3, {1,0,0,
			       0,1,0,
			       0,0,1}),
		generate_addmod_worms(3),
	};*/
}	


inline int site_basis::size() const {
	return states.size();
}
