#pragma once

#include <cstdint>
#include <cassert>
#include "bond.h"

using jm_action = uint8_t;

class jm {
public:
	static const int basis_size = 4;
	static int local_basis_size(const site &s);

	uint8_t code() const;

	double j(const site &s) const;
	double m(const site &s) const;

	char name(const site &s) const;
	
	jm_action to(jm target) const;
	jm apply(jm_action action);

	jm() = default;
	explicit jm(uint8_t code);
	bool operator==(const jm &other) const;
private:
	uint8_t code_{};
};

inline jm::jm(uint8_t code)
	: code_{code} {
}

inline char jm::name(const site &s) const {
	if(s.nspinhalfs == 1) {
		switch(code_) {
		case 0:
			return '+';
		case 1:
			return '-';
		default:
			return 'x';
		}
	} else if(s.nspinhalfs == 2) {
		switch(code_) {
		case 0:
			return '*';
		case 1:
			return 'P';
		case 2:
			return '0';
		case 3:
			return 'M';
		default:
			assert(false);
		}
	}

	return 'X';
}


inline bool jm::operator==(const jm &other) const {
	return code_ == other.code_;
}

inline uint8_t jm::code() const {
	return code_;
}

inline jm_action jm::to(jm target) const {
	return code_ ^ target.code_;
}

inline jm jm::apply(jm_action action) {
	code_ ^= action;
	return *this;
}

inline int jm::local_basis_size(const site &s) {
	if(s.nspinhalfs == 1) {
		return 2;
	}
	return 4;
}

inline double jm::j(const site &s) const {
	if(s.nspinhalfs == 1) {
		return 0.5;
	} else if(s.nspinhalfs == 2) {
		if(code_ == 0) {
			return 0;
		}
		return 1;
	}
	assert(false);
}

inline double jm::m(const site &s) const {
	if(s.nspinhalfs == 1) {
		switch(code_) {
		case 0:
			return 0.5;
		case 1: 
			return -0.5;
		default:
			assert(false);
		}
	} else if(s.nspinhalfs == 2) {
		switch(code_) {
		case 0:
			return 0;
		case 1: 
			return 1;
		case 2:
			return 0;
		case 3:
			return -1;
		default:
			assert(false);
		}
	}
	assert(false);
}

