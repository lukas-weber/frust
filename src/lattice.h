#pragma once

#include "bond.h"
#include "opercode.h"
#include "vertex_data.h"
#include <nlohmann/json.hpp>
#include <vector>

struct unitcell {
	vec2 a1;
	vec2 a2;

	std::vector<uc_site> sites;
	std::vector<uc_bond> bonds;
};

class lattice {
public:
	unitcell uc;
	std::vector<int> uc_signs;

	int Lx{}, Ly{};
	int spinhalf_count{};
	double energy_offset{};

	std::vector<lat_site> sites;
	std::vector<lat_bond> bonds;

	const vertex_data &get_vertex_data(int bond) const;
	const uc_site &get_uc_site(int site) const;

	auto split_idx(int idx) const;

	void vertex_print() const;

	lattice(const unitcell &uc, const std::vector<int> &uc_signs, int Lx, int Ly, bool with_vertex_data);

	void to_json(nlohmann::json &out);

private:
	std::vector<vertex_data> vertices_; // [uc_bond]

	void calculate_energy_offset();
	void init_vertex_data(const unitcell &uc);
};

inline auto lattice::split_idx(int idx) const {
	int iuc = idx % uc.sites.size();
	idx /= uc.sites.size();

	int x = idx % Lx;
	int y = idx / Lx;

	return std::tuple{iuc, x, y};
}

inline const vertex_data &lattice::get_vertex_data(int bond) const {
	assert(vertices_.size());
	return vertices_[bond % vertices_.size()];
}
inline const uc_site &lattice::get_uc_site(int site) const {
	return uc.sites[site % uc.sites.size()];
}
