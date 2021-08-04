#pragma once

#include <array>
#include <vector>
#include "basis.h"

using vec2 = std::array<double,2>;

inline const vec2 operator+(const vec2 &a, const vec2 &b) {
	vec2 r;
	for(int i = 0; i < 2; i++) {
		r[i] = a[i] + b[i];
	}
	return r;
}

inline const vec2 operator*(double a, const vec2 &b) {
	vec2 r;

	for(int i = 0; i < 2; i++) {
		r[i] = a * b[i];
	}
	return r;
}

struct uc_bond {
	int i{};
	struct {
		int dx{};
		int dy{};
		int uc{};
	} j;
	
	// bonds between individual constituent spins
	// shape should be [spin_i * site_j.nspinhalfs + spin_j]
	std::vector<double> J;
};

struct uc_site {
	vec2 pos;

	// bonds inside the side. vector enumerates the ascendingly ordered pairs (i, j) with i < j.
	std::vector<double> Jin;
	site_basis basis;

	// magnetic field
	double h{};

	int coordination{}; // filled in automatically
};

struct lat_bond {
	int i;
	int j;
};

struct lat_site {
	vec2 pos;
};
