#pragma once
#include <cstdint>
#include <fmt/format.h>

#include "bond.h"
#include "jm.h"

class opercode {
public:
	using vertex_leg = uint32_t;
	using vertex_idx = uint32_t;

	static opercode make_identity();
	static opercode make_vertex(uint32_t bond, jm s0, jm s1, jm s2, jm s3);

	std::string name(const site &si, const site &sj) const;

	uint32_t code() const;
	vertex_idx vertex() const;
	int bond() const;

	jm_action action(vertex_leg l) const;
	jm leg_state(vertex_leg l) const;

	bool identity() const;
	bool diagonal() const;

	opercode() = default;
	explicit opercode(uint32_t code)
		: code_{code} {
	}
private:
	uint32_t code_{};
};

inline std::string opercode::name(const site &si, const site &sj) const {
	if(identity()) {
		return "[ident]";
	}
	return fmt::format("[{}{}>{}{}]", leg_state(0).name(si), leg_state(1).name(sj), leg_state(2).name(si), leg_state(3).name(sj));
}


inline jm_action opercode::action(vertex_leg side) const {
	assert(side == 0 || side == 1);
	return (vertex()>>(2*side+4))&0x3;
}

inline jm opercode::leg_state(vertex_leg leg) const {
	assert(leg < 4);
	if(leg == 0 || leg == 1) {
		return jm{static_cast<uint8_t>((vertex()>>(2*leg))&0x3)};
	}

	return leg_state(leg&1).apply(action(leg&1));
}

inline opercode opercode::make_identity() {
	return opercode{0};
}

inline opercode opercode::make_vertex(uint32_t bond, jm s0, jm s1, jm s2, jm s3) {
	jm_action s02 = s0.to(s2);
	jm_action s13 = s1.to(s3);

	assert(bond < (1<<19));
	
	return opercode{1 | s0.code() << 1 | s1.code() << 3 | s02 << 5 | s13 << 7 | bond << 9};
}

inline uint32_t opercode::code() const {
	return code_;
}

inline int opercode::bond() const {
	return code_>>9;
}

inline opercode::vertex_idx opercode::vertex() const {
	return (code_&((1<<9)-1))>>1;
}

inline bool opercode::identity() const {
	return code_ == 0;
}

inline bool opercode::diagonal() const {
	return (vertex() >> 4) == 0;
}

