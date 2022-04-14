#pragma once

#include "util/vec.h"
#include <nlohmann/json.hpp>

struct unitcell {
	struct bond {
		int i{};
		struct {
			int dx{};
			int dy{};
			int uc{};
		} j;
	};

	struct site {
		vec2 pos;
		int sublattice_sign{1};
		int coordination{}; // filled automatically
	};

	vec2 a1;
	vec2 a2;

	std::vector<site> sites;
	std::vector<bond> bonds;
};

struct lattice {
	struct bond {
		int type{};
		int i{};
		int j{};
	};
	unitcell uc;

	int Lx{}, Ly{};

	std::vector<bond> bonds;

	auto split_idx(int site_idx) const;
	vec2 site_pos(int site_idx) const;
	double site_sublattice_sign(int site_idx) const;
	int site_count() const;
	void to_json(nlohmann::json &out) const;

	lattice(const unitcell &uc, int Lx, int Ly);
};

inline auto lattice::split_idx(int site_idx) const {
	int iuc = site_idx % uc.sites.size();
	site_idx /= uc.sites.size();

	int x = site_idx % Lx;
	int y = site_idx / Lx;

	return std::tuple{iuc, x, y};
}

inline int lattice::site_count() const {
	return Lx * Ly * uc.sites.size();
}

inline vec2 lattice::site_pos(int site_idx) const {
	auto [x, y, iuc] = split_idx(site_idx);
	return (x + uc.sites[iuc].pos[0]) * uc.a1 + (y + uc.sites[iuc].pos[1]) * uc.a2;
}

inline double lattice::site_sublattice_sign(int site_idx) const {
	return uc.sites[site_idx % uc.sites.size()].sublattice_sign;
}

namespace unitcells {
const unitcell square{{1, 0}, {0, 1}, {{{0, 0}}}, {{0, {1, 0, 0}}, {0, {0, 1, 0}}}};

const unitcell columnar_dimer{{1, 0},
                              {0, 2},
                              {{{0, 0}, 1}, {{0, 0.5}, -1}},
                              {{0, {0, 0, 1}}, {0, {1, 0, 0}}, {1, {1, 0, 1}}, {1, {0, 1, 0}}}};

const unitcell shastry_sutherland{{1, 0},
                                  {0, 1},
                                  {{{0, 0}}, {{0.5, 0}}, {{0.5, 0.5}}, {{0, 0.5}}},
                                  {{0, {0, 0, 1}},
                                   {1, {0, 0, 2}},
                                   {2, {0, 0, 3}},
                                   {3, {0, 0, 0}},
                                   {1, {0, 0, 3}},
                                   {1, {1, 0, 0}},
                                   {2, {1, 0, 3}},
                                   {2, {0, 1, 1}},
                                   {3, {0, 1, 0}},
                                   {2, {1, 1, 0}}}};

const unitcell shastry_sutherland_dimer{
    {1, 0},
    {0, 1},
    {{{0, 0}}, {{0.5, 0.5}}},
    {{0, {0, 0, 1}}, {1, {1, 0, 0}}, {1, {0, 1, 0}}, {1, {1, 1, 0}}}};

const unitcell triangle{
    {1, 0}, {0, 1}, {{{0, 0}}}, {{0, {0, 1, 0}}, {0, {1, 1, 0}}, {0, {1, 0, 0}}}};

const unitcell kagome{{1, 0},
                      {-0.5, sqrt(3) / 2},
                      {{{0, 0}}, {{0.5, 0}}, {{0.5, 0.5}}},
                      {{0, {0, 0, 1}},
                       {1, {0, 0, 2}},
                       {2, {0, 0, 0}},
                       {1, {1, 0, 0}},
                       {2, {1, 1, 0}},
                       {2, {0, 1, 1}}}};

const unitcell kagome_dimer_spin{
    kagome.a1,
    kagome.a2,
    {{{0, 0}}, {{0.2, 0.5}}},
    {{0, {0, 0, 1}}, {0, {1, 0, 0}}, {1, {0, 1, 0}}, {1, {1, 1, 0}}},
};

const unitcell kagome_trimer = triangle;

const unitcell lieb_lattice{{1, 0},
                            {0, 1},
                            {{{0, 0}, 1}, {{0.5, 0}, -1}, {{0, 0.5}, -1}},
                            {{0, {0, 0, 1}}, {0, {0, 0, 2}}, {1, {1, 0, 0}}, {2, {0, 1, 0}}}};

const unitcell triangle_square{{1, 0},
                               {0, 2},
                               {{{-0.3, 0}}, {{0.3, 0}}, {{0, 0.5}}},
                               {{0, {0, 0, 1}},
                                {1, {0, 0, 2}},
                                {2, {0, 0, 0}},
                                {1, {1, 0, 0}},
                                {2, {0, 1, 0}},
                                {2, {0, 1, 1}}}};

const unitcell triangle_square_dimer_spin{triangle_square.a1,
                                          triangle_square.a2,
                                          {{{0, 0}}, {{0, 0.5}}},
                                          {{0, {0, 0, 1}}, {0, {1, 0, 0}}, {1, {0, 1, 0}}}};

const unitcell triangle_square_dimer_spin_rot{
    triangle_square.a1,
    triangle_square.a2,
    {{{0.5, 0.5}}, {{0, 0}}},
    {{0, {0, 0, 1}}, {0, {1, 0, 1}}, {0, {0, 1, 0}}, {0, {0, 1, 1}}}};

const unitcell triangle_square_trimer = square;
}
