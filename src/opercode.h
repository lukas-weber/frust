#pragma once
#include <cassert>
#include <cstdint>
#include <fmt/format.h>

using opercode_uint = uint64_t;

class vertexcode {
public:
	static const opercode_uint maxbits = 8 * 3 + 1;
	bool diagonal() const {
		return code_ & 1L;
	}

	vertexcode() = default;

	explicit vertexcode(bool diagonal, opercode_uint vertex_idx)
	    : code_{diagonal | (vertex_idx << 1)} {
		assert(!invalid());
	}

	explicit vertexcode(opercode_uint code) : code_{code} {}

	opercode_uint vertex_idx() const {
		return code_ >> 1L;
	}

	opercode_uint code() const {
		return code_;
	}

	bool invalid() const {
		return code_ > 1L << maxbits;
	}

private:
	opercode_uint code_{(1L << maxbits) + 1L};
};

class opercode {
public:
	static opercode make_identity();

	opercode_uint code() const;
	vertexcode vertex() const;
	opercode_uint bond() const;

	bool identity() const;
	bool diagonal() const;

	opercode() = default;
	explicit opercode(opercode_uint code) : code_{code} {}
	explicit opercode(opercode_uint bond, vertexcode vertex);

private:
	opercode_uint code_{};
};

inline opercode opercode::make_identity() {
	return opercode{0};
}

inline opercode::opercode(opercode_uint bond, vertexcode vertex) {
	opercode_uint v = vertex.code();

	assert(v < (1L << vertexcode::maxbits));
	assert(bond < (1L << (8L * sizeof(code_) - vertexcode::maxbits - 1L)));

	code_ = 1L | v << 1L | bond << (1L + vertexcode::maxbits);
}

inline opercode_uint opercode::code() const {
	return code_;
}

inline opercode_uint opercode::bond() const {
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
