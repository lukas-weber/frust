#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "../src/downfolded_peierls_coupling.cpp"

using namespace downfolded_peierls_coupling;

TEST_CASE("coupling evaluation") {
	// REQUIRE(disp_op(0, 4, 0.2).real() == Approx(0.0003201315461539436));
	// REQUIRE(disp_op(2, 4, 0.2).real() == Approx(-0.0661082733373851));
	// REQUIRE(disp_op(2, 3, 0.2).imag() == Approx(0.3260592963812132));
	// REQUIRE(disp_op(1, 3, 0.01).real() == Approx(-0.00012246428128910691));

	CHECK(generator({{0.55, 0.5, 10}}).elem(1, 1) == Approx(1.0951812739150684));
	CHECK(generator({{0.23, 0.001, 10}}).elem(3, 1) == Approx(1.9189577301773e-7));

	CHECK(generator({{0.55, 0.5, 10}}).elem(1, 1) == Approx(1.0951812739150684));
	CHECK(generator({{0.23, 0.2, 10}}).elem(3, 3) == Approx(1.005934546456162));
	CHECK(generator({{0.55, 0.5, 10}}).elem(0, 3) == 0);
	CHECK(generator({{0.23, 0.2, 10}}).elem(1, 3) == Approx(0.007657606067842092));
	CHECK(generator({{0.23, 0.2, 10}}).elem(3, 1) ==
	      Approx(generator({{0.23, 0.2, 4}}).elem(3, 1)));

	CHECK(generator({{0.54, 0.14, 10}, {0.54, 0.14, 10}}).elem(33, 11) ==
	      Approx(-0.0198192171961849));
	CHECK(generator({{0.54, 0.14, 10}, {0.54, 0.14, 10}}).elem(33, 15) ==
	      Approx(-0.015138320606728406));
}

TEST_CASE("stability") {
	for(int n = 1; n < 100; n += 10) {
		CHECK_NOTHROW(generator({{0.49, 1, n}}).matrix());
	}
}
