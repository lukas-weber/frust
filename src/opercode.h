#pragma once
#include <cassert>
#include <cstdint>
#include <fmt/format.h>

class vertexcode {
public:
	static const uint64_t maxbits = 8 * 3 + 1;
	bool diagonal() const {
		return code_ & 1L;
	}

	vertexcode() = default;

	explicit vertexcode(bool diagonal, uint64_t vertex_idx) : code_{diagonal | (vertex_idx << 1)} {
		assert(!invalid());
	}

	explicit vertexcode(uint64_t code) : code_{code} {}

	uint64_t vertex_idx() const {
		return code_ >> 1L;
	}

	uint64_t code() const {
		return code_;
	}

	bool invalid() const {
		return code_ > 1L << maxbits;
	}

private:
	uint64_t code_{(1L << maxbits) + 1L};
};

class opercode {
public:
	static opercode make_identity();

	uint64_t code() const;
	vertexcode vertex() const;
	uint64_t bond() const;

	bool identity() const;
	bool diagonal() const;

	opercode() = default;
	explicit opercode(uint64_t code) : code_{code} {}
	explicit opercode(uint64_t bond, vertexcode vertex);

private:
	uint64_t code_{};
};

inline opercode opercode::make_identity() {
	return opercode{0};
}

inline opercode::opercode(uint64_t bond, vertexcode vertex) {
	uint64_t v = vertex.code();

	assert(v < (1L << vertexcode::maxbits));
	assert(bond < (1L << (8L * sizeof(code_) - vertexcode::maxbits - 1L)));

	code_ = 1L | v << 1L | bond << (1L + vertexcode::maxbits);
}

inline uint64_t opercode::code() const {
	return code_;
}

inline uint64_t opercode::bond() const {
	return code_ >> (1L + vertexcode::maxbits);
}

inline vertexcode opercode::vertex() const {
	return vertexcode{(code_ & ((1L << vertexcode::maxbits) - 1L)) >> 1L};
}

inline bool opercode::identity() const {
	return code_ == 0L;
}

inline bool opercode::diagonal() const {
	return vertex().diagonal();
}
