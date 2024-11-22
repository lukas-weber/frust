#include <catch2/catch.hpp>

#include "models/cavity_magnet/downfolded_coupling.cpp"

TEST_CASE("coupling evaluation") {
	REQUIRE(dispOp(0, 4, 0.2).real() == Approx(0.0003201315461539436));
	REQUIRE(dispOp(2, 4, 0.2).real() == Approx(-0.0661082733373851));
	REQUIRE(dispOp(2, 3, 0.2).imag() == Approx(0.3260592963812132));
	REQUIRE(dispOp(1, 3, 0.01).real() == Approx(-0.00012246428128910691));

	CHECK(J(1, 1, {{0.55, 0.5, 10}}) == Approx(1.0951812739150684));
	CHECK(J(3, 1, {{0.23, 0.001, 10}}) == Approx(1.9189577301773e-7));

	CHECK(J(1, 1, {{0.55, 0.5, 10}}) == Approx(1.0951812739150684));
	CHECK(J(3, 3, {{0.23, 0.2, 10}}) == Approx(1.005934546456162));
	CHECK(J(0, 3, {{0.55, 0.5, 10}}) == 0);
	CHECK(J(1, 3, {{0.23, 0.2, 10}}) == Approx(0.007657606067842092));
	CHECK(J(3, 1, {{0.23, 0.2, 10}}) == Approx(J(3, 1, {{0.23, 0.2, 4}})));

	CHECK(J(33, 11, {{0.54, 0.14, 10}, {0.54, 0.14, 10}}) == Approx(-0.0198192171961849));
	CHECK(J(33, 15, {{0.54, 0.14, 10}, {0.54, 0.14, 10}}) == Approx(-0.015138320606728406));
}
