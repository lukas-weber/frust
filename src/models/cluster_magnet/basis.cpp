#include "basis.h"
#include "../common/spinop.h"

site_basis::site_basis(int nspinhalfs, const std::vector<state> &states,
                       const Eigen::MatrixXd &trans)
    : states_{states}, nspinhalfs{nspinhalfs}, trans_{trans} {}

static const std::array<Eigen::MatrixXd, 2> Shalf = spin_operators(2);

static Eigen::MatrixXd lift_spinop(const Eigen::Matrix2d &op, int nspinhalfs, int pos) {
	assert(pos < nspinhalfs);
	int dimleft = 1 << pos;
	int dimright = 1 << (nspinhalfs - pos - 1);

	return kronecker_prod(kronecker_prod(Eigen::MatrixXd::Identity(dimleft, dimleft), op),
	                      Eigen::MatrixXd::Identity(dimright, dimright));
}

std::array<Eigen::MatrixXd, 2> site_basis::spinop(int spinhalf) const {
	assert(spinhalf < nspinhalfs);
	std::array<Eigen::MatrixXd, 2> transformed_spin;
	std::transform(Shalf.begin(), Shalf.end(), transformed_spin.begin(), [&](const auto &op) {
		auto lifted = lift_spinop(op, nspinhalfs, spinhalf);
		Eigen::MatrixXd result = trans_.transpose() * lifted * trans_;

		return result;
	});

	return transformed_spin;
}
