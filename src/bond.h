#pragma once

#include "basis.h"
#include <Eigen/Dense>

using vec2 = Eigen::Vector2d;

struct uc_bond {
	int i{};
	struct {
		int dx{};
		int dy{};
		int uc{};
	} j;
};

struct uc_site {
	vec2 pos;
	int coordination{}; // filled automatically
};

// FIXME: remove these?
struct lat_bond {
	int i;
	int j;
};

struct lat_site {
	vec2 pos;
	int dim;
};
