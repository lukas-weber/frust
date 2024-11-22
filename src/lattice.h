#pragma once

#include <vector>
#include <nlohmann/json.hpp>
#include "bond.h"
#include "vertex_data.h"
#include "opercode.h"

struct unitcell {
	vec2 a1;
	vec2 a2;

	std::vector<uc_site> sites;
	std::vector<uc_bond> bonds;
};


class lattice {
public:
	unitcell uc;

	int Lx{}, Ly{};
	int spinhalf_count{};
	double energy_offset{};

	std::vector<lat_site> sites;
	std::vector<lat_bond> bonds;

	const vertex_data &get_vertex_data(int bond) const;
	const uc_site &get_uc_site(int site) const;
	opercode vertex_idx_opercode(int bond, int vertex_idx) const;

	void vertex_print() const;

	lattice(const unitcell &uc, int Lx, int Ly);

	void to_json(nlohmann::json &out);
private:
	std::vector<vertex_data> vertices_; // [uc_bond]

	void calculate_energy_offset();
	void init_sublattice();
	void init_vertex_data(const unitcell &uc);
};

inline const vertex_data &lattice::get_vertex_data(int bond) const {
	return vertices_[bond % vertices_.size()];
}
inline const uc_site &lattice::get_uc_site(int site) const {
	return uc.sites[site%uc.sites.size()];
}

inline opercode lattice::vertex_idx_opercode(int bond, int vertex_idx) const {
	const auto &ls = vertices_[bond%vertices_.size()].legstates[vertex_idx];

	return opercode::make_vertex(bond, ls[0], ls[1], ls[2], ls[3]);
}
