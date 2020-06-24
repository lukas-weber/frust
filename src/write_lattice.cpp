#include "write_lattice.h"
#include "latticedef.h"
#include <loadleveller/loadleveller.h>
#include <nlohmann/json.hpp>

void write_lattice(int argc, char **argv) {
	if(argc != 3) {
		printf("Error: Invalid number of arguments.\nUsage: isingsse lattice JOBFILE TASKNAME\n\n");
		printf(
		    "This exports lattice data in a human readable format containing all bond and spin "
		    "properties.\n");
		return;
	}

	std::string jobfilename{argv[1]};
	std::string taskname{argv[2]};

	loadl::parser jobfile{jobfilename};
	auto task = jobfile["tasks"][taskname];

	lattice lat = lattice_from_param(task, true);

	nlohmann::json out;
	lat.to_json(out);
	std::cout << out.dump(1) << std::endl;
}

