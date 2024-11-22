#include "frust.h"
#include "loadleveller/loadleveller.h"
//#include "write_lattice.h"
#include "vertices.h"

int main(int argc, char *argv[]) {
	if(argc > 1 && std::string(argv[1]) == "lattice") {
//		write_lattice(argc - 1, argv + 1);
		return 0;
	}
	vertex_data v{bond{0, 1, 1.2, 1, 1}, site{{0,0}, 2}, site{{0,0},2}};
	v.print();
	
	//loadl::run<frust>(argc, argv);
}
