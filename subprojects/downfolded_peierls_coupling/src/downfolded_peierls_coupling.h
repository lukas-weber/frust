#pragma once
#include <complex>
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

class generator {
public:
	generator(const std::vector<mode_params> &modes, double tolerance = 1e-10);

	double elem(int n, int m) const;
	std::vector<double> matrix() const;

private:
	std::complex<double> disp_op(int n, int l, double g) const;

	double tolerance_{};
	double logtolerance_{};
	int lmax_{};

	std::vector<mode_params> modes_;
	std::vector<double> lgamma_cache_;
};

}
