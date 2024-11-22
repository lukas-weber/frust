#include "downfolded_peierls_coupling.h"
#include "occupation_numbers.h"
#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <complex>
#include <numeric>

namespace downfolded_peierls_coupling {

using cmplx = std::complex<double>;

generator::generator(const std::vector<mode_params> &modes, double tolerance)
    : tolerance_{tolerance}, modes_(modes) {
	int max_max_photons = std::transform_reduce(
	    modes.begin(), modes.end(), 0, [](int a, int b) { return std::max(a, b); },
	    [](const mode_params &m) { return m.max_photons; });

	lmax_ = 2 * std::max(static_cast<int>(modes.size()), max_max_photons) - log(tolerance);

	lgamma_cache_.resize(lmax_ + 1);

	for(int i = 0; i < static_cast<int>(lgamma_cache_.size()); i++) {
		lgamma_cache_[i] = lgamma(1 + i);
	}
}

cmplx generator::disp_op(int n, int l, double g) const {
	using namespace std::complex_literals;
	if(g == 0) {
		return n == l;
	}

	cmplx sum = 0;
	double logg = log(g);

	std::array<cmplx, 4> sign = {1, 1i, -1, -1i};

	for(int a = 0; a <= std::min(n, l); a++) {
		int pow_exp = n + l - 2 * a;
		double term = exp(logg * pow_exp + 0.5 * lgamma_cache_[l] + 0.5 * lgamma_cache_[n] -
		                  lgamma_cache_[a] - lgamma_cache_[n - a] - lgamma_cache_[l - a]);
		sum += sign[pow_exp % sign.size()] * term;
	}
	return sum;
}

double generator::elem(int n, int m) const {
	double sum = 0;
	double term = 0;

	std::vector<int> mode_dims(modes_.size());
	std::transform(modes_.begin(), modes_.end(), mode_dims.begin(),
	               [](const auto &m) { return m.max_photons; });
	std::vector<int> l_dims(modes_.size(), lmax_);
	occupation_numbers nit{n, mode_dims}, mit{m, mode_dims};

	for(int l = 0; l < round(std::pow(lmax_, modes_.size())); l++) {
		cmplx disp_prod = 1;
		double denom_n = 1;
		double denom_m = 1;
		occupation_numbers lit{l, l_dims};
		int mode_idx{};
		for(auto ni = nit.begin(), mi = mit.begin(), li = lit.begin();
		    ni != nit.end() && mi != mit.end() && li != lit.end(); ++ni, ++mi, ++li) {
			const auto &mo = modes_[mode_idx];

			disp_prod *= disp_op(*ni, *li, mo.g) * std::conj(disp_op(*mi, *li, mo.g));
			denom_n += mo.omega * (*li - *ni);
			denom_m += mo.omega * (*li - *mi);
			mode_idx++;
		}
		term = disp_prod.real() * (1.0 / denom_n + 1.0 / denom_m);
		sum += term;
	}
	assert(fabs(term) < tolerance_);
	int sign = std::inner_product(nit.begin(), nit.end(), mit.begin(), 1, std::multiplies(),
	                              [](int ni, int mi) { return 1 - 2 * ((ni - mi) / 2 & 1); });
	double g2_sum = std::transform_reduce(modes_.begin(), modes_.end(), 0.0, std::plus(),
	                                      [](const auto &m) { return m.g * m.g; });

	return exp(-g2_sum) * sign * 0.5 * sum;
}

std::vector<double> generator::matrix() const {
	int photon_dim = std::transform_reduce(modes_.begin(), modes_.end(), 1, std::multiplies(),
	                                       [](const auto &m) { return m.max_photons; });

	std::vector<double> coupling(photon_dim * photon_dim);
	for(int m = 0; m < photon_dim; m++) {
		for(int n = 0; n < photon_dim; n++) {
			coupling[m * photon_dim + n] = elem(m, n);
		}
	}
	return coupling;
}

}
