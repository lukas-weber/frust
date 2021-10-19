#include "frust.h"
#include "loadleveller/loadleveller.h"
#include "write_lattice.h"

int main(int argc, char *argv[]) {
	if(argc > 1 && std::string(argv[1]) == "lattice") {
		write_lattice(argc - 1, argv + 1);
		return 0;
	}

	loadl::run<frust>(argc, argv);
}
