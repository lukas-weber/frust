#pragma once
#include <cstdint>
#include <fmt/format.h>
#include <cassert>

class vertexcode {
public:
	static const uint32_t maxbits = 4 * 3 + 1;
	bool diagonal() const {
		return code_ & 1;
	}

	vertexcode() = default;

	explicit vertexcode(bool diagonal, uint32_t vertex_idx) : code_{diagonal | (vertex_idx << 1)} {}

	explicit vertexcode(uint32_t code) : code_{code} {}

	uint32_t vertex_idx() const {
		return code_ >> 1;
	}

	uint32_t code() const {
		return code_;
	}

	bool invalid() const {
		return code_ > 1 << maxbits;
	}

private:
	uint32_t code_{(1 << maxbits) + 1};
};

class opercode {
public:
	static opercode make_identity();

	uint32_t code() const;
	vertexcode vertex() const;
	uint32_t bond() const;

	bool identity() const;
	bool diagonal() const;

	opercode() = default;
	explicit opercode(uint32_t code) : code_{code} {}
	explicit opercode(uint32_t bond, vertexcode vertex);

private:
	uint32_t code_{};
};

inline opercode opercode::make_identity() {
	return opercode{0};
}

inline opercode::opercode(uint32_t bond, vertexcode vertex) {
	uint32_t v = vertex.code();

	assert(v < (1 << vertexcode::maxbits));
	assert(bond < (1 << (8 * sizeof(code_) - vertexcode::maxbits - 1)));

	code_ = 1 | v << 1 | bond << (1 + vertexcode::maxbits);
}

inline uint32_t opercode::code() const {
	return code_;
}

inline uint32_t opercode::bond() const {
	return code_ >> (1 + vertexcode::maxbits);
}

inline vertexcode opercode::vertex() const {
	return vertexcode{(code_ & ((1 << (1 + 4 * 3)) - 1)) >> 1};
}

inline bool opercode::identity() const {
	return code_ == 0;
}

inline bool opercode::diagonal() const {
	return vertex().diagonal();
}
