#pragma once
#include <cstdint>
#include <fmt/format.h>

#include "bond.h"
#include "basis.h"

class opercode {
public:
	using vertex_leg = uint32_t;
	using vertex_idx = uint32_t;

	static opercode make_identity();
	static opercode make_vertex(uint32_t bond, state_idx s0, state_idx s1, state_idx s2, state_idx s3);

	std::string name(const site_basis &si, const site_basis &sj) const;

	uint32_t code() const;
	vertex_idx vertex() const;
	int bond() const;

	state_idx leg_state(vertex_leg l) const;

	bool identity() const;
	bool diagonal() const;

	opercode() = default;
	explicit opercode(uint32_t code)
		: code_{code} {
	}
private:
	uint32_t code_{};
};

inline std::string opercode::name(const site_basis &bi, const site_basis &bj) const {
	if(identity()) {
		return "[ident]";
	}
	return fmt::format("[{}{}>{}{}]", bi.states[leg_state(0)].name, bj.states[leg_state(1)].name, bi.states[leg_state(2)].name, bj.states[leg_state(3)].name);
}

inline state_idx opercode::leg_state(vertex_leg leg) const {
	assert(leg < 4);
	return (vertex()>>(site_basis::state_bits*leg))&((1<<site_basis::state_bits)-1);
}

inline opercode opercode::make_identity() {
	return opercode{0};
}

inline opercode opercode::make_vertex(uint32_t bond, state_idx s0, state_idx s1, state_idx s2, state_idx s3) {
	const uint32_t bits = site_basis::state_bits;
	assert(bond < (1<<(8*sizeof(code_)-1-4*bits)));
	return opercode{1 | s0 << 1 | s1 << (1+bits) | s2 << (1+2*bits) | s3 << (1+3*bits) | bond << (1+4*bits)};
}

inline uint32_t opercode::code() const {
	return code_;
}

inline int opercode::bond() const {
	return code_>>(1+4*site_basis::state_bits);
}

inline opercode::vertex_idx opercode::vertex() const {
	return (code_&((1<<(1+4*site_basis::state_bits))-1))>>1;
}

inline bool opercode::identity() const {
	return code_ == 0;
}

inline bool opercode::diagonal() const {
	return leg_state(0) == leg_state(2) && leg_state(1) == leg_state(3);
}

