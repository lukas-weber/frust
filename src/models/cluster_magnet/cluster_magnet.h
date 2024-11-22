#pragma once

#include "../common/lattice.h"
#include "../model.h"
#include "basis.h"
#include "measurement_settings.h"
#include <nlohmann/json.hpp>
#include <optional>

struct cluster_bond {
	// bonds between individual constituent spins
	// shape should be [spin_i * site_j.nspinhalfs + spin_j]
	std::vector<double> J;
};

struct cluster_site {
	// bonds inside the cluster site. vector enumerates the ascendingly ordered pairs (i, j) with i
	// < j.
	std::vector<double> Jin;
	const site_basis &basis;

	// magnetic field
	double h{};
};

class cluster_magnet : public model {
public:
	lattice lat;
	cluster_magnet_measurement_settings settings;
	int spinhalf_count{};

	cluster_magnet(const lattice &lat, const std::vector<cluster_site> &sites,
	               const std::vector<cluster_bond> &bonds,
	               const cluster_magnet_measurement_settings &settings);

	const cluster_bond &get_bond(int bond_idx) const;
	const cluster_site &get_site(int site_idx) const;

	const site_basis &get_basis(int site_idx) const;
	std::optional<lattice::site_idx> get_lattice_site_idx(int site_idx) const;

	void register_evalables(loadl::evaluator &eval, double T) const override;

	sse_data generate_sse_data() const override;
	void to_json(nlohmann::json &out) const override;
	int normalization_site_count() const override;

private:
	std::vector<cluster_bond> bonds_;
	std::vector<cluster_site> sites_;
};

inline const cluster_bond &cluster_magnet::get_bond(int bond_idx) const {
	return bonds_[bond_idx % bonds_.size()];
}

inline const cluster_site &cluster_magnet::get_site(int site_idx) const {
	return sites_[site_idx % sites_.size()];
}

inline const site_basis &cluster_magnet::get_basis(int site_idx) const {
	return sites_[site_idx % sites_.size()].basis;
}

inline std::optional<lattice::site_idx> cluster_magnet::get_lattice_site_idx(int site_idx) const {
	return site_idx;
}

inline int cluster_magnet::normalization_site_count() const {
	return spinhalf_count;
}
