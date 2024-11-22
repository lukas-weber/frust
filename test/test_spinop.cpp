#include "models/common/spinop.h"
#include "util/kronecker_product.h"
#include <catch2/catch.hpp>

TEST_CASE("spinop generation") {
	auto shalfops = spin_operators(2);

	// clang-format off
	auto shalfz = init_mat(2, 2, {0.5,0,0,-0.5,});
	auto shalfplus = init_mat(2, 2, {0,1,0,0,});
	// clang-format on

	REQUIRE(shalfops[0] == shalfplus);
	REQUIRE(shalfops[1] == shalfz);
}
