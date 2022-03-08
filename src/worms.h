#pragma once
#include <cstdint>

// Implementation of the worms. Choice is largely irrelevant as long as it is bijective in a certain sense.
using worm_idx = int;
using state_idx = uint_fast8_t;
/*
inline state_idx worm_action(worm_idx worm, state_idx state, int basis_size) {
	return (state + worm + 1) % basis_size;
}

inline worm_idx worm_inverse(worm_idx worm, int basis_size) {
	return basis_size - worm - 2;
}

inline int worm_count(int basis_size) {
	return basis_size - 1;
}*/
inline state_idx worm_action(worm_idx worm, state_idx state, int) {
	return state^(worm + 1);
}

inline worm_idx worm_inverse(worm_idx worm, int) {
	return worm;
}

inline int worm_count(int basis_size) {
	return basis_size - 1;
}
	
