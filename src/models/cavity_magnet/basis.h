#pragma once

#include "occupation_numbers.h"
#include "worms.h" // XXX: move state_idx type to another header?
#include <functional>
#include <numeric>

class cavity_basis {
public:
	enum class site_type {
		spin,
		mode,
	};

	double n(state_idx state, int mode) const {
		if(type_ == site_type::spin) {
			return 0;
		}

		return photon_numbers_[state * mode_num_ + mode];
	}

	double m(state_idx state) const {
		if(type_ == site_type::mode) {
			return 0;
		}
		return (spin_dim_ - 1) * 0.5 - state;
	}

	static cavity_basis make_spin_basis(int spin_dim) {
		return cavity_basis{site_type::spin, spin_dim, 0, {}};
	}

	static cavity_basis make_mode_basis(const std::vector<int> &mode_dims) {
		int photon_dim = std::accumulate(mode_dims.begin(), mode_dims.end(), 1, std::multiplies());
		std::vector<int> photon_numbers;
		int mode_num = mode_dims.size();
		photon_numbers.reserve(photon_dim * mode_num);
		for(int n = 0; n < photon_dim; n++) {
			for(int d : occupation_numbers{n, mode_dims}) {
				photon_numbers.push_back(d);
			}
		}

		return cavity_basis{site_type::mode, 0, mode_num, std::move(photon_numbers)};
	}

private:
	cavity_basis(site_type type, int spin_dim, int mode_num, std::vector<int> &&photon_numbers)
	    : type_{type}, spin_dim_{spin_dim}, mode_num_{mode_num}, photon_numbers_{photon_numbers} {}

	site_type type_;
	int spin_dim_{};
	int mode_num_{};
	std::vector<int> photon_numbers_; // [state_idx; mode_num]
};
