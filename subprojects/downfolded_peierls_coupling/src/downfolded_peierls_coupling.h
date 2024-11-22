#pragma once
#include <vector>

namespace downfolded_peierls_coupling {

struct mode_params {
	double omega{};
	double g{};
	int max_photons{};

	// implementation of constructor is needed for the python bindings
	mode_params(double omega, double g, int max_photons)
	    : omega{omega}, g{g}, max_photons{max_photons} {}
	mode_params() = default;
};

double elem(int n, int m, const std::vector<mode_params> &modes, double tolerance);
std::vector<double> matrix(const std::vector<mode_params> &modes);

};
