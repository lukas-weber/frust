#pragma once

#include "basis.h"

// Implementation of the worms. Choice is largely irrelevant as long as it is bijective in a certain sense.
using worm_idx = int;

inline state_idx worm_action(worm_idx worm, state_idx state, int basis_size) {
	return (state + worm + 1) % basis_size;
}

inline worm_idx worm_inverse(worm_idx worm, int basis_size) {
	return basis_size - worm - 2;
}

inline int worm_count(int basis_size) {
	return basis_size - 1;
}
	
