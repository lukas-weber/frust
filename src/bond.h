#pragma once

#include <array>
#include <vector>

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

struct bond {
	int i{0};
	int j{0};
	
	// bonds between individual constituent spins
	// shape should be [spin_i * site_j.nspinhalfs + spin_j]
	std::vector<double> J;
};

struct site {
	vec2 pos;

	int nspinhalfs{}; // S=1/2 spins contained in this site
	double Jin{};
	
	int8_t sublattice{};
};
