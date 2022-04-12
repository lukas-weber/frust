#include "downfolded_coupling.h"
#include <cmath>
#include <complex>

using cmplx = std::complex<double>;
using namespace std::complex_literals;

static cmplx dispOp(int n, int l, double g) {
	cmplx sum = 0;

	for(int a = 0; a <= std::min(n, l); a++) {
		cmplx term = std::pow(1i * g, n + l - 2 * a) *
		             exp(0.5 * lgamma(1 + l) + 0.5 * lgamma(1 + n) - lgamma(1 + a) -
		                 lgamma(1 + n - a) - lgamma(1 + l - a));
		sum += term;
	}
	return sum;
}

double J(int n, int m, double omega, double g, double tolerance = 1e-12) {
	if((n + m) % 2 != 0) {
		return 0;
	}
	double sum = 0;
	double term = 0;
	for(int l = 0; l < std::max(n, m) + 2 || fabs(term) > tolerance; l++) {
		term = (dispOp(n, l, g) * std::conj(dispOp(m, l, g))).real() *
		       (1.0 / (1.0 + omega * (l - n)) + 1.0 / (1.0 + omega * (l - m)));
		sum += term;
	}

	return (1 - 2 * ((n - m) / 2 & 1)) * 0.5 * exp(-g * g) * sum;
}

Eigen::MatrixXd downfolded_coupling(double omega, double g, int max_bosons) {
	Eigen::MatrixXd coupling(max_bosons, max_bosons);
	for(int m = 0; m < max_bosons; m++) {
		for(int n = 0; n < max_bosons; n++) {
			coupling(m, n) = J(n, m, omega, g);
		}
	}
	return coupling;
}
