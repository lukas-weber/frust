#pragma once
#include <Eigen/Dense>

struct downfolded_coupling_params {
	double omega{};
	double g{};
	int max_photons{};
};

Eigen::MatrixXd downfolded_coupling(const std::vector<downfolded_coupling_params> &modes);
