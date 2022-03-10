#include "spinop.h"

std::array<Eigen::MatrixXd, 2> spin_operators(int dimension) {
	double S = (dimension - 1) / 2.0;
	Eigen::MatrixXd sz = Eigen::MatrixXd::Zero(dimension, dimension);
	Eigen::MatrixXd splus = Eigen::MatrixXd::Zero(dimension, dimension);

	for(int i = 0; i < dimension; i++) {
		double m = S - i;
		sz(i, i) = m;
		if(i < dimension - 1) {
			splus(i, i + 1) = sqrt((S - m + 1) * (S + m));
		}
	}

	return {splus, sz};
}
