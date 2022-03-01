#pragma once

#include "model.h"
#include "util/lattice.h"
#include <nlohmann/json.hpp>

struct cluster_bond {
	// bonds between individual constituent spins
	// shape should be [spin_i * site_j.nspinhalfs + spin_j]
	std::vector<double> J;
};

struct cluster_site {
	// bonds inside the cluster site. vector enumerates the ascendingly ordered pairs (i, j) with i < j.
	std::vector<double> Jin;
	site_basis basis;

	// magnetic field
	double h{};
	int sublattice_sign{1};
};


class cluster_magnet : public model {
public:
	lattice lat;
	int spinhalf_count{};

	cluster_magnet(const lattice &lat, const std::vector<cluster_site>& sites, const std::vector<cluster_bond>& bonds);

	const cluster_bond& get_bond(int bond_idx) const;
	const cluster_site& get_site(int site_idx) const;

	sse_data generate_sse_data() const override;
	void to_json(nlohmann::json& out) const override;
	int normalization_site_count() const override;
private:
	std::vector<cluster_bond> bonds_;
	std::vector<cluster_site> sites_;
};

inline const cluster_bond& cluster_magnet::get_bond(int bond_idx) const {
	return bonds_[bond_idx % bonds_.size()];
}

inline const cluster_site& cluster_magnet::get_site(int site_idx) const {
	return sites_[site_idx % sites_.size()];
}

inline int cluster_magnet::normalization_site_count() const {
	return spinhalf_count;
}
