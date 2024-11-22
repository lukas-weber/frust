#include <catch2/catch.hpp>

#include "models/cavity_magnet/downfolded_coupling.cpp"

TEST_CASE("coupling evaluation") {
	REQUIRE(dispOp(0, 4, 0.2).real() == Approx(0.0003265986323710905));
	REQUIRE(dispOp(2, 4, 0.2).real() == Approx(-0.06744374904565532));
	REQUIRE(dispOp(2, 3, 0.2).imag() == Approx(0.3326461310962948));

	CHECK(J(1, 1, 0.55, 0.5) == Approx(1.0951812739150684));
	CHECK(J(3, 3, 0.23, 0.2) == Approx(1.005934546456162));
	CHECK(J(0, 3, 0.55, 0.5) == 0);
	CHECK(J(1, 3, 0.23, 0.2) == Approx(-0.007657606067842092));
	CHECK(J(3, 1, 0.23, 0.2) == Approx(J(3, 1, 0.23, 0.2)));
}
