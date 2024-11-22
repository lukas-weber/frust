#include "model_def.h"
#include "cluster_magnet_def.h"
#include "cavity_magnet_def.h"

std::unique_ptr<model> model_from_param(const loadl::parser &p) {
	std::string model_name = p.get<std::string>("model", "cluster_magnet"); // TODO: remove default when backwards compatibility is broken.
	
	if(model_name == "cluster_magnet") {
		return cluster_magnet_from_param(p);
	}
	if(model_name == "cavity_magnet") {
		return cavity_magnet_from_param(p);
	}

	throw std::runtime_error{fmt::format("unknown model '{}'", model_name)};
}
