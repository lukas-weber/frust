#pragma once
#include <vector>

struct downfolded_coupling_params {
	double omega{};
	double g{};
	int max_photons{};
};

std::vector<double> downfolded_coupling(const std::vector<downfolded_coupling_params> &modes);
