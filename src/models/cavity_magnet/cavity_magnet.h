#pragma once

#include "../common/lattice.h"
#include "../model.h"
#include "basis.h"
#include <optional>

class cavity_magnet : public model {
public:
	struct site {
		int spin_dim{2}; // = 2 * S + 1
	};

	struct bond {
		double J{1};
	};

	struct mode {
		double omega{};
		double coupling{};
		int max_bosons{};
	};

	lattice lat;

	cavity_magnet(const lattice &lat, const std::vector<mode> &modes,
	              const std::vector<site> &sites, const std::vector<bond> &bonds);

	const bond &get_bond(int bond_idx) const {
		return bonds_[bond_idx % bonds_.size()];
	}

	cavity_basis get_basis(int abstract_site_idx) const {
		return bases_[abstract_site_idx];
	}

	std::optional<lattice::site_idx> get_lattice_site_idx(int abstract_site_idx) const {
		if(abstract_site_idx < static_cast<int>(modes_.size())) {
			return std::nullopt;
		}
		return abstract_site_idx - modes_.size();
	}

	sse_data generate_sse_data() const override;
	void to_json(nlohmann::json &out) const override;
	int normalization_site_count() const override {
		return lat.site_count();
	}

private:
	std::vector<mode> modes_;
	std::vector<site> sites_;
	std::vector<bond> bonds_;
	std::vector<cavity_basis> bases_;
};
