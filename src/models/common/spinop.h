#pragma once

#include <Eigen/Dense>
#include <array>

// spin_operators returns {splus, sz} for an arbitrary spin magnitude specified by dimension = 2*S+1
std::array<Eigen::MatrixXd, 2> spin_operators(int dimension);
