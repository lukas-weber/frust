#include "vertices.h"

#include <Eigen/Dense>

template<int d>
static Eigen::Matrix<double, d * d, d * d> kronecker_prod(const Eigen::Matrix<double, d, d> &a,
                                                          const Eigen::Matrix<double, d, d> &b) {
	Eigen::Matrix<double, d * d, d * d> res;
	for(int i1 = 0; i1 < d; i1++) {
		for(int j1 = 0; j1 < d; j1++) {
			for(int i2 = 0; i2 < d; i2++) {
				for(int j2 = 0; j2 < d; j2++) {
					int i = d * i1 + i2;
					int j = d * j1 + j2;
					res(i, j) = a(i1, j1) * b(i2, j2);
				}
			}
		}
	}

	return res;
}

vertex_data::vertex_data(const bond &b) {
	using mat = Eigen::Matrix<double, 6, 6>;

	std::vector<mat> S(3);

	S[0].diagonal(1) << 0, 1, 0, sqrt(2), sqrt(2);
	S[1] = S[0].transpose();
	S[2].diagonal() << 0, 0.5, -0.5, 1, 0, -1;

	auto H = b.J * (0.5 * (kronecker_prod(S[0], S[1]) + kronecker_prod(S[1], S[0])) +
	                kronecker_prod(S[2], S[2]));
	mat S2 = 0.5 * (S[0] * S[1] + S[1] * S[0]) + S[2] * S[2];
	mat Id = mat::Identity();
	H += b.Ki * kronecker_prod(S2, Id);
	H += b.Kj * kronecker_prod(Id, S2);
}
