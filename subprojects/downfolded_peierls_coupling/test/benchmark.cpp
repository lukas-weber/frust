#include "downfolded_peierls_coupling.h"

int main(int, char **) {
	downfolded_peierls_coupling::generator gen({{0.49, 0.99, 100}});
	gen.matrix();
}
