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
	void apply(jm_action action);

	jm() = default;
	explicit jm(uint32_t code);
	bool operator==(const jm &other) const;
private:
	uint32_t code_{};
};

inline jm::jm(uint32_t code)
	: code_{code} {
}

inline bool jm::operator==(const jm &other) const {
	return code_ == other.code_;
}

inline uint32_t jm::code() const {
	return code_;
}

inline jm_action jm::to(jm target) const {
	return code_ ^ target.code_;
}

inline void jm::apply(jm_action action) {
	code_ ^= action;
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

