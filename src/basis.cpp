#include "basis.h"
#include "util.h"

static std::array<Eigen::Matrix2d, 3> gen_spinhalfop() {
	std::array<Eigen::Matrix2d, 3> S;
	S[0] = Eigen::Matrix2d::Zero();
	S[1] = Eigen::Matrix2d::Zero();
	S[2] = Eigen::Matrix2d::Zero();
	
	S[0].diagonal(1) << 1;
	S[1] = S[0].transpose();
	S[2].diagonal() << 0.5, -0.5;
	return S;
}


static const std::array<Eigen::Matrix2d, 3> Shalf = gen_spinhalfop();

static Eigen::MatrixXd lift_spinop(const Eigen::Matrix2d &op, int nspinhalfs, int pos) {
	assert(pos < nspinhalfs);
	int dimleft = 1<<pos;
	int dimright = 1<<(nspinhalfs-pos-1);

	return kronecker_prod(kronecker_prod(Eigen::MatrixXd::Identity(dimleft,dimleft), op), Eigen::MatrixXd::Identity(dimright,dimright));
}

std::array<Eigen::MatrixXd,3> site_basis::spinop(int spinhalf) const {

	assert(spinhalf < nspinhalfs);
	std::array<Eigen::MatrixXd,3> transformed_spin;
	std::transform(Shalf.begin(), Shalf.end(), transformed_spin.begin(), [&](const auto &op) {
		auto lifted = lift_spinop(op, nspinhalfs, spinhalf);
		Eigen::MatrixXd result = trans_.transpose() * lifted * trans_;

		return result; 
	});

	return transformed_spin;
}
	