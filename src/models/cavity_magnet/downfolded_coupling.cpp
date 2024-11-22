#include "downfolded_coupling.h"
#include <cmath>
#include <complex>
#include <numeric>

using cmplx = std::complex<double>;
using namespace std::complex_literals;

static cmplx dispOp(int n, int l, double g) {
	if(g == 0) {
		return n == l;
	}

	cmplx sum = 0;

	for(int a = 0; a <= std::min(n, l); a++) {
		cmplx term = std::pow(1i * g, n + l - 2 * a) *
		             exp(0.5 * lgamma(1 + l) + 0.5 * lgamma(1 + n) - lgamma(1 + a) -
		                 lgamma(1 + n - a) - lgamma(1 + l - a));
		sum += term;
	}
	return exp(-g * g / 2) * sum;
}

double J(int n, int m, const std::vector<downfolded_coupling_params> &modes,
         double tolerance = 1e-12) {
	double sum = 0;
	double term = 0;
	int lmax = -log(tolerance); // convergence should be exponential, at least for g < 1

	for(int l = 0; l < round(std::pow(lmax, modes.size())); l++) {
		cmplx disp_prod = 1;
		double denom_n = 1;
		double denom_m = 1;
		int dim = 1;
		int ldim = 1;
		for(const auto &mo : modes) {
			int ni = (n / dim) % mo.max_photons;
			int mi = (m / dim) % mo.max_photons;
			int li = (l / ldim) % lmax;
			dim *= mo.max_photons;
			ldim *= lmax;

			disp_prod *= dispOp(ni, li, mo.g) * std::conj(dispOp(mi, li, mo.g));
			denom_n += mo.omega * (li - ni);
			denom_m += mo.omega * (li - mi);
		}
		term = disp_prod.real() * (1.0 / denom_n + 1.0 / denom_m);
		sum += term;
	}
	assert(fabs(term) < tolerance);
	int sign = 1;
	int dim = 1;
	for(const auto &mo : modes) {
		int ni = (n / dim) % mo.max_photons;
		int mi = (m / dim) % mo.max_photons;
		dim *= mo.max_photons;

		sign *= (1 - 2 * ((ni - mi) / 2 & 1));
	}
	return sign * 0.5 * sum;
}

Eigen::MatrixXd downfolded_coupling(const std::vector<downfolded_coupling_params> &modes) {
	int photon_dim = std::accumulate(modes.begin(), modes.end(), 1,
	                                 [](int a, const auto &m) { return a * m.max_photons; });

	Eigen::MatrixXd coupling(photon_dim, photon_dim);
	for(int m = 0; m < photon_dim; m++) {
		for(int n = 0; n < photon_dim; n++) {
			coupling(m, n) = J(m, n, modes);
		}
	}
	return coupling;
}
