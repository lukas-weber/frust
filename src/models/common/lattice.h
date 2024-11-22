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
	using site_idx = int;

	struct bond {
		int type{};
		int i{};
		int j{};
	};
	unitcell uc;

	int Lx{}, Ly{};

	std::vector<bond> bonds;

	auto split_idx(site_idx) const;
	vec2 site_pos(site_idx) const;
	double site_sublattice_sign(site_idx) const;
	int site_count() const;
	void to_json(nlohmann::json &out) const;

	lattice(const unitcell &uc, int Lx, int Ly);
};

inline auto lattice::split_idx(site_idx i) const {
	int iuc = i % uc.sites.size();
	i /= uc.sites.size();

	int x = i % Lx;
	int y = i / Lx;

	return std::tuple{iuc, x, y};
}

inline int lattice::site_count() const {
	return Lx * Ly * uc.sites.size();
}

inline vec2 lattice::site_pos(site_idx i) const {
	auto [iuc, x, y] = split_idx(i);
	Eigen::Matrix2d A;
	A.col(0) = uc.a1;
	A.col(1) = uc.a2;
	return A * (vec2{x, y} + uc.sites[iuc].pos);
}

inline double lattice::site_sublattice_sign(site_idx i) const {
	return uc.sites[i % uc.sites.size()].sublattice_sign;
}

namespace unitcells {
const unitcell square{{1, 0}, {0, 1}, {{{0, 0}}}, {{0, {1, 0, 0}}, {0, {0, 1, 0}}}};

const unitcell columnar_dimer{{1, 0},
                              {0, 2},
                              {{{0, 0}, 1}, {{0, 0.5}, -1}},
                              {{0, {0, 0, 1}}, {0, {1, 0, 0}}, {1, {1, 0, 1}}, {1, {0, 1, 0}}}};

const unitcell honeycomb{{sqrt(3) / 2, -0.5},
                         {sqrt(3) / 2, 0.5},
                         {{{0, 0}, 1}, {{1.0 / 3, 1.0 / 3}, -1}},
                         {{0, {0, 0, 1}}, {1, {0, 1, 0}}, {1, {1, 0, 0}}}};

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
