#pragma once

#include "../common/lattice.h"
#include "../model.h"

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

	sse_data generate_sse_data() const override;
	void to_json(nlohmann::json &out) const override;
	int normalization_site_count() const override {
		return lat.site_count();
	}

private:
	std::vector<mode> modes_;
	std::vector<site> sites_;
	std::vector<bond> bonds_;
};
