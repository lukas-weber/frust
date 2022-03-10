#pragma once

#include <Eigen/Dense>
#include <array>

template<typename MatType1, typename MatType2>
Eigen::MatrixXd kronecker_prod(const MatType1 &a, const MatType2 &b) {
	Eigen::MatrixXd res(a.rows() * b.rows(), a.cols() * b.cols());
	for(int i1 = 0; i1 < a.rows(); i1++) {
		for(int j1 = 0; j1 < a.cols(); j1++) {
			for(int i2 = 0; i2 < b.rows(); i2++) {
				for(int j2 = 0; j2 < b.cols(); j2++) {
					int i = b.rows() * i1 + i2;
					int j = b.cols() * j1 + j2;
					res(i, j) = a(i1, j1) * b(i2, j2);
				}
			}
		}
	}

	return res;
}

// computes the “scalar product” for two sets of spin operators {Splus, Sz}
template<typename MatrixIn>
Eigen::MatrixXd scalar_product(const std::array<MatrixIn, 2> &Si,
                               const std::array<MatrixIn, 2> &Sj) {
	return 0.5 * (kronecker_prod(Si[0], Sj[0].transpose()) +
	              kronecker_prod(Si[0].transpose(), Sj[0])) +
	       kronecker_prod(Si[1], Sj[1]);
}

inline Eigen::MatrixXd init_mat(int rows, int cols, std::initializer_list<double> elems) {
	int idx{};
	Eigen::MatrixXd res(rows, cols);
	for(auto e : elems) {
		int i = idx / cols;
		int j = idx % cols;
		res(i, j) = e;
		idx++;
	}
	return res;
}
