#pragma once
#include <cstdint>

#include "bond.h"
#include "jm.h"

class opercode {
public:
	using vertex_leg = uint32_t;

	static const uint32_t max_vertex_count = 64;

	static opercode make_identity();
	static opercode make_vertex(uint32_t bond, jm s0, jm s1, jm s2, jm s3);

	uint32_t code() const;
	int vertex() const;
	int bond() const;

	jm_action action(vertex_leg l) const;

	bool identity() const;
	bool diagonal() const;

	opercode() = default;
	explicit opercode(uint32_t code)
		: code_{code} {
	}
private:
	uint32_t code_{};
};


inline jm_action opercode::action(vertex_leg side) const {
	return (vertex()-1)>>(3*side+6);
}

inline opercode opercode::make_identity() {
	return opercode{0};
}

inline opercode opercode::make_vertex(uint32_t bond, jm s0, jm s1, jm s2, jm s3) {
	jm_action s02 = s0.to(s2);
	jm_action s13 = s1.to(s3);

	assert(bond < (1<<19));
	
	return opercode{1 | s0.code() << 1 | s1.code() << 4 | s02 << 7 | s13 << 10 | bond << 13};
}

inline uint32_t opercode::code() const {
	return code_;
}

inline int opercode::bond() const {
	return code_>>13;
}

inline int opercode::vertex() const {
	return code_&((1<<11)-1);
}

inline bool opercode::identity() const {
	return code_ == 0;
}

inline bool opercode::diagonal() const {
	return vertex() >> 7 == 0;
}

