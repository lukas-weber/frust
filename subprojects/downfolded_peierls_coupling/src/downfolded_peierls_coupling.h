#pragma once
#include <vector>

struct downfolded_coupling_params {
	double omega{};
	double g{};
	int max_photons{};

	// implementation of constructor is needed for the python bindings
	downfolded_coupling_params(double omega, double g, int max_photons)
	    : omega{omega}, g{g}, max_photons{max_photons} {}
	downfolded_coupling_params() = default;
};

std::vector<double> downfolded_coupling(const std::vector<downfolded_coupling_params> &modes);
