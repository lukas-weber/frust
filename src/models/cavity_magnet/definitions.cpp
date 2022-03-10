#include "definitions.h"
#include "cavity_magnet.h"

struct cavity_magnet_proto {
	unitcell uc;
	std::vector<cavity_magnet::site> sites;
	std::vector<cavity_magnet::bond> bonds;
};

static cavity_magnet_proto make_square(const loadl::parser &p) {
	double J = p.get<double>("J");
	int spin_dim = round(2 * p.get<double>("S", 0.5) + 1);

	return {unitcells::square, {{spin_dim}}, {{J}, {J}}};
}

std::unique_ptr<cavity_magnet> cavity_magnet_from_param(const loadl::parser &p) {
	auto lat = p.get<std::string>("lattice");

	cavity_magnet_proto proto;
	int Lx = p.get<int>("Lx");
	int Ly = p.get<int>("Ly", Lx);

	if(lat == "square") {
		proto = make_square(p);
	}

	int max_bosons = p.get<int>("max_bosons");

	auto freqs = p.get<std::vector<double>>("mode_freqs");
	auto couplings = p.get<std::vector<double>>("mode_couplings");

	if(freqs.size() != couplings.size()) {
		throw std::runtime_error{
		    fmt::format("mode_freqs ({}) and mode_couplings ({}) do not match in length",
		                freqs.size(), couplings.size())};
	}

	std::vector<cavity_magnet::mode> modes;
	for(int i = 0; i < static_cast<int>(freqs.size()); i++) {
		modes.push_back({freqs[i], couplings[i], max_bosons});
	}

	return std::make_unique<cavity_magnet>(lattice{proto.uc, Lx, Ly}, modes, proto.sites,
	                                       proto.bonds);
}
