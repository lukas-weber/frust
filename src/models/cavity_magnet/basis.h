#pragma once

#include "worms.h" // XXX: move state_idx type to another header?

class cavity_basis {
public:
	enum class site_type {
		spin,
		mode,
	};

	double n(state_idx state) const {
		if(type_ == site_type::spin) {
			return 0;
		}
		return state;
	}

	double m(state_idx state) const {
		if(type_ == site_type::mode) {
			return 0;
		}
		return (spin_dim_ - 1) * 0.5 - state;
	}

	static cavity_basis make_spin_basis(int spin_dim) {
		return cavity_basis{site_type::spin, spin_dim};
	}

	static cavity_basis make_mode_basis() {
		return cavity_basis{site_type::mode, 0};
	}

private:
	cavity_basis(site_type type, int spin_dim) : type_{type}, spin_dim_{spin_dim} {}
	site_type type_;
	int spin_dim_{};
};
