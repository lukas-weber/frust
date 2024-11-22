#pragma once

#include <cstdint>
#include <cassert>

using jm_action = uint32_t;

class jm {
public:
	uint32_t code() const;

	double j() const;
	double m() const;
	
	jm_action to(jm target) const;
	jm apply(jm_action action) const;

	explicit jm(uint32_t code);
private:
	uint32_t code_{};
};

jm::jm(uint32_t code)
	: code_{code} {
}

inline uint32_t jm::code() const {
	return code_;
}

inline jm_action jm::to(jm target) const {
	return code_ ^ target.code_;
}

inline jm jm::apply(jm_action action) const {
	return jm{code_ ^ action};
}

inline double jm::j() const {
	switch(code_) {
	case 0:
		return 0;
	case 1:
		return 0.5;
	case 2:
		return 0.5;
	case 3:
		return 1;
	case 4:
		return 1;
	case 5:
		return 1;
	default:
		assert(false);
		return 0;
	}
}

inline double jm::m() const {
	switch(code_) {
	case 0:
		return 0;
	case 1:
		return 0.5;
	case 2:
		return -0.5;
	case 3:
		return 1;
	case 4:
		return 0;
	case 5:
		return -1;
	default:
		assert(false);
		return 0;
	}
}

