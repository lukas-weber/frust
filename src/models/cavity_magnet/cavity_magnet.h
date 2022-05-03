#pragma once

#include "../common/lattice.h"
#include "../model.h"
#include "basis.h"
#include "measurement_settings.h"
#include <optional>

class cavity_magnet : public model {
public:
	struct site {
		int spin_dim{2}; // = 2 * S + 1
	};

	struct bond {
		double J{1};
		std::vector<double> mode_couplings;
	};

	struct mode {
		double omega{};
		int max_photons{};
	};

	lattice lat;
	std::vector<mode> modes;
	cavity_magnet_measurement_settings settings;

	cavity_magnet(const lattice &lat, const std::vector<mode> &modes,
	              const std::vector<site> &sites, const std::vector<bond> &bonds, double U,
	              const cavity_magnet_measurement_settings &settings);

	const bond &get_bond(int bond_idx) const {
		return bonds_[bond_idx % bonds_.size()];
	}

	cavity_basis get_basis(int abstract_site_idx) const {
		return bases_[abstract_site_idx];
	}

	std::optional<lattice::site_idx> get_lattice_site_idx(int abstract_site_idx) const {
		if(abstract_site_idx == 0) {
			return std::nullopt;
		}
		return abstract_site_idx - 1;
	}

	void register_evalables(loadl::evaluator &eval, double T) const override;

	sse_data generate_sse_data() const override;
	void to_json(nlohmann::json &out) const override;
	int normalization_site_count() const override {
		return lat.site_count();
	}

private:
	std::vector<site> sites_;
	std::vector<bond> bonds_;
	std::vector<cavity_basis> bases_;

	double U_{}; // Hubbard interaction
};
