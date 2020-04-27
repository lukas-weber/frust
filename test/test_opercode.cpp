#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include "opercode.h"

TEST_CASE("opercode") {
	int bond = GENERATE(take(10, random(0, 100 * 100)));

	uint8_t basis_size = 1<<site_basis::state_bits;

	for(uint8_t i0 = 0; i0 < basis_size; i0++) {
		for(uint8_t i1 = 0; i1 < basis_size; i1++) {
			for(uint8_t i2 = 0; i2 < basis_size; i2++) {
				for(uint8_t i3 = 0; i3 < basis_size; i3++) {
					state_idx l0{i0};
					state_idx l1{i1};
					state_idx l2{i2};
					state_idx l3{i3};

					opercode op = opercode::make_vertex(bond, l0, l1, l2, l3);

					CHECK(op.leg_state(0) == l0);
					CHECK(op.leg_state(1) == l1);
					CHECK(op.leg_state(2) == l2);
					CHECK(op.leg_state(3) == l3);
					CHECK(op.bond() == bond);

					CHECK(op.diagonal() == (l0 == l2 && l1 == l3));
					CHECK(!op.identity());
				}
			}
		}
	}

	CHECK(opercode::make_identity().identity());
}
