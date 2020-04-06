#pragma once

#include <vector>
#include <nlohmann/json.hpp>
#include "bond.h"
#include "vertex_data.h"
#include "opercode.h"

struct unitcell {
	vec2 a1;
	vec2 a2;

	std::vector<site> sites;
	std::vector<bond> bonds;
};


class lattice {
public:
	int Lx{}, Ly{};
	int spinhalf_count{};
	double energy_offset{};

	std::vector<site> sites;
	std::vector<bond> bonds;

	const vertex_data::transition &vertex_transition(opercode op, int leg_in, jm_action action_in) const;
	double vertex_weight(opercode op) const;
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

inline const vertex_data::transition &lattice::vertex_transition(opercode op, int leg_in, jm_action action_in) const {
	return vertices_[op.bond()%vertices_.size()].get_transition(op, leg_in, action_in);
}

inline double lattice::vertex_weight(opercode op) const {
	return vertices_[op.bond()%vertices_.size()].get_weight(op);
}

inline opercode lattice::vertex_idx_opercode(int bond, int vertex_idx) const {
	const auto &ls = vertices_[bond%vertices_.size()].legstates[vertex_idx];

	return opercode::make_vertex(bond, ls[0], ls[1], ls[2], ls[3]);
}
