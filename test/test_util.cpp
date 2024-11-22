#include "util/kronecker_product.h"
#include <catch2/catch.hpp>

TEST_CASE("init_mat") {
	Eigen::MatrixXd a(3, 3);
	a << 0, 1, 2, 3, 4, 5, 6, 7, 8;

	REQUIRE(a == init_mat(3, 3, {0, 1, 2, 3, 4, 5, 6, 7, 8}));
}
