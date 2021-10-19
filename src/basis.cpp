#include "basis.h"
#include "util.h"
#include <fmt/format.h>

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

std::vector<worm_function> generate_xor_worms(uint32_t basis_size) {
	std::vector<worm_function> res;
	for(uint32_t i = 1; i < basis_size; i++) {
		std::vector<state_idx> action;
		for(uint32_t j = 0; j < basis_size; j++) {
			action.push_back(j ^ i);
		}
		res.push_back(worm_function{"^", static_cast<int>(i - 1), action});
	}
	return res;
}

std::vector<worm_function> generate_addmod_worms(uint32_t basis_size) {
	std::vector<worm_function> res;
	for(uint32_t i = 1; i < basis_size; i++) {
		std::vector<state_idx> action;
		for(uint32_t j = 0; j < basis_size; j++) {
			action.push_back((j + i) % basis_size);
		}
		res.push_back(worm_function{"+", static_cast<int>(basis_size - i - 1), action});
	}
	return res;
}

site_basis::site_basis(int nspinhalfs, const std::vector<state> &states,
                       const Eigen::MatrixXd &trans, const std::vector<worm_function> &worms)
    : nspinhalfs{nspinhalfs}, states{states}, worms{worms}, trans_{trans} {}

static const std::array<Eigen::Matrix2d, 3> Shalf = gen_spinhalfop();

static Eigen::MatrixXd lift_spinop(const Eigen::Matrix2d &op, int nspinhalfs, int pos) {
	assert(pos < nspinhalfs);
	int dimleft = 1 << pos;
	int dimright = 1 << (nspinhalfs - pos - 1);

	return kronecker_prod(kronecker_prod(Eigen::MatrixXd::Identity(dimleft, dimleft), op),
	                      Eigen::MatrixXd::Identity(dimright, dimright));
}

std::array<Eigen::MatrixXd, 3> site_basis::spinop(int spinhalf) const {
	assert(spinhalf < nspinhalfs);
	std::array<Eigen::MatrixXd, 3> transformed_spin;
	std::transform(Shalf.begin(), Shalf.end(), transformed_spin.begin(), [&](const auto &op) {
		auto lifted = lift_spinop(op, nspinhalfs, spinhalf);
		Eigen::MatrixXd result = trans_.transpose() * lifted * trans_;

		return result;
	});

	return transformed_spin;
}
