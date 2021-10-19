#define CATCH_CONFIG_MAIN
#include "opercode.h"
#include <catch2/catch.hpp>

TEST_CASE("opercode") {
	uint32_t bond = GENERATE(take(10, random(0, 100 * 100)));

	for(uint32_t vidx = 0; vidx < 77; vidx++) {
		for(int diagonal = 0; diagonal < 2; diagonal++) {
			vertexcode v{!!diagonal, vidx};
			CHECK(v.diagonal() == diagonal);
			CHECK(v.vertex_idx() == vidx);

			opercode op{bond, v};

			CHECK(op.bond() == bond);
			CHECK(op.diagonal() == diagonal);
			CHECK(op.vertex().code() == v.code());

			CHECK(!op.identity());
		}
	}

	CHECK(opercode::make_identity().identity());
}
