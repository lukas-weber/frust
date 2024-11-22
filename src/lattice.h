#pragma once

#include <vector>
#include <nlohmann/json.hpp>
#include "bond.h"

struct unitcell {
	vec2 a1;
	vec2 a2;

	std::vector<site> sites;
	std::vector<bond> bonds;
};


class lattice {
public:
	int Lx{}, Ly{};

	std::vector<int8_t> sublattice;

	std::vector<site> sites;
	std::vector<bond> bonds;

	double energy_offset() const;

	lattice(unitcell uc, int Lx, int Ly);

	void to_json(nlohmann::json &out);
private:
	double energy_offset_{};

	void calculate_energy_offset();
	void init_sublattice();
};
