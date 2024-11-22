#pragma once

#include "bond.h"
#include <vector>

class vertex_data {
public: 
	static const int basis_size = 6;
	using leg = uint8_t;
	
	vertex_data(const bond &b);
private:
	std::vector<std::array<leg,4>> legstates_;
	std::vector<double> weights_;
};
