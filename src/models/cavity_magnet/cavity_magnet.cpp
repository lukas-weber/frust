#include "cavity_magnet.h"
#include "../common/mag_est.h"
#include "../common/spinop.h"
#include "downfolded_coupling.h"
#include "nlohmann/json.hpp"
#include "photon_est.h"
#include "util/kronecker_product.h"
#include <functional>

cavity_magnet::cavity_magnet(const lattice &lat, const std::vector<mode> &modes,
                             const std::vector<site> &sites, const std::vector<bond> &bonds,
                             const cavity_magnet_measurement_settings &settings)
    : model{model_type::cavity_magnet}, lat{lat}, modes{modes}, settings{settings}, sites_{sites},
      bonds_{bonds} {
	assert(bonds_.size() == lat.uc.bonds.size());
	assert(sites_.size() == lat.uc.sites.size());
	for(const auto &m : modes) {
		(void)m;
		bases_.push_back(cavity_basis::make_mode_basis());
	}

	for(int i = 0; i < lat.site_count(); i++) {
		bases_.push_back(cavity_basis::make_spin_basis(sites[i % sites.size()].spin_dim));
	}
}

sse_data cavity_magnet::generate_sse_data() const {
	std::vector<vertex_data> vert_data;
	std::vector<sse_data::bond> sse_bonds;
	std::vector<sse_data::site> sse_sites;

	int photon_dim = std::transform_reduce(modes.begin(), modes.end(), 1, std::multiplies{},
	                                       [](const auto &mode) { return mode.max_photons; });

	Eigen::MatrixXd Hphot = Eigen::MatrixXd::Zero(photon_dim, photon_dim);
	for(int n = 0; n < photon_dim; n++) {
		int dim = 1;
		for(const auto &m : modes) {
			int ni = (n / dim) % m.max_photons;
			Hphot(n, n) += m.omega * ni;
		}
	}

	int mode_count = modes.size();

	int bond_idx{};
	for(const auto &b : bonds_) {
		int mode_idx{};

		std::vector<downfolded_coupling_params> mode_params(modes.size());
		for(const auto &m : modes) {
			mode_params[mode_idx] = {m.omega, b.mode_couplings[mode_idx], m.max_photons};
			mode_idx++;
		}

		int site_idx_i = lat.uc.bonds[bond_idx].i;
		int site_idx_j = lat.uc.bonds[bond_idx].j.uc;

		int spin_dim_i = sites_[site_idx_i].spin_dim;
		int spin_dim_j = sites_[site_idx_j].spin_dim;

		auto spinop_i = spin_operators(spin_dim_i);
		auto spinop_j = spin_operators(spin_dim_j);

		Eigen::MatrixXd spin_identity =
		    Eigen::MatrixXd::Identity(spin_dim_i * spin_dim_j, spin_dim_i * spin_dim_j);

		Eigen::MatrixXd exchange_photon_coupling = downfolded_coupling(mode_params);

		Eigen::MatrixXd H =
		    b.J * kronecker_prod(exchange_photon_coupling,
		                         -0.25 * spin_identity + scalar_product(spinop_i, spinop_j)) +
		    1.0 / lat.bonds.size() * kronecker_prod(Hphot, spin_identity);

		vert_data.push_back({{photon_dim, spin_dim_i, spin_dim_j}, H});
		bond_idx++;
	}

	bond_idx = 0;
	for(const auto &b : lat.bonds) {
		for(int mode_idx = 0; mode_idx < static_cast<int>(modes.size()); mode_idx++) {
			int bond_type = b.type * mode_count + mode_idx;
			sse_bonds.push_back({bond_type, {mode_idx, mode_count + b.i, mode_count + b.j}});
		}
		bond_idx++;
	}

	std::transform(modes.begin(), modes.end(), std::back_inserter(sse_sites),
	               [](const auto &m) { return sse_data::site{m.max_photons}; });
	for(int i = 0; i < lat.Lx * lat.Ly; i++) {
		for(int uc = 0; uc < static_cast<int>(lat.uc.sites.size()); uc++) {
			sse_sites.push_back({sites_[uc].spin_dim});
		}
	}

	return sse_data{vert_data, sse_sites, sse_bonds};
}

void cavity_magnet::to_json(nlohmann::json &out) const {
	out["model"] = "cavity_magnet";
	lat.to_json(out);

	int bond_idx{};
	for(const auto &b : lat.bonds) {
		(void)b;
		out["bonds"][bond_idx]["J"] = get_bond(bond_idx).J;
		out["bonds"][bond_idx]["mode_couplings"] = get_bond(bond_idx).mode_couplings;
		bond_idx++;
	}

	int site_idx{};
	for(int i = 0; i < lat.Lx * lat.Ly; i++) {
		for(int uc = 0; uc < static_cast<int>(lat.uc.sites.size()); uc++) {
			out["sites"][site_idx]["spin_dim"] = sites_[uc].spin_dim;
			site_idx++;
		}
	}

	for(const auto &m : modes) {
		out["modes"].push_back({{"omega", m.omega}, {"max_photons", m.max_photons}});
	}
}

void cavity_magnet::register_evalables(loadl::evaluator &eval, double T) const {
	using M = cavity_magnet;
	if(settings.measure_mag) {
		mag_est<mag_sign::none, M>{*this, T, 0}.register_evalables(eval);
	}

	if(settings.measure_sxsymag) {
		mag_est<mag_sign::x | mag_sign::y, M>{*this, T, 0}.register_evalables(eval);
	}

	if(settings.measure_sxsymag) {
		mag_est<mag_sign::x, M>{*this, T, 0}.register_evalables(eval);
	}

	if(settings.measure_sxsucmag) {
		mag_est<mag_sign::x | mag_sign::uc, M>{*this, T, 0}.register_evalables(eval);
	}

	photon_est{*this, 0}.register_evalables(eval);
}
