#include "definitions.h"
#include "cavity_magnet.h"
#include "measurement_settings.h"

enum class cavity_type {
	polarized,  // omega, g, phi: single, linearly polarized mode
	degenerate, // omega, g: two degenerate modes with perpendicular polarization
};

cavity_type cavity_from_name(const std::string &name) {
	if(name == "polarized") {
		return cavity_type::polarized;
	}
	if(name == "degenerate") {
		return cavity_type::degenerate;
	}
	throw std::runtime_error{fmt::format("unknown cavity type: {}", name)};
}

struct cavity_magnet_proto {
	unitcell uc;
	std::vector<cavity_magnet::site> sites;
	std::vector<cavity_magnet::bond> bonds;
};

static cavity_magnet_proto make_square(const loadl::parser &p) {
	double J = p.get<double>("J");
	int spin_dim = round(2 * p.get<double>("S", 0.5) + 1);

	std::vector<double> mode_coupling_x;
	std::vector<double> mode_coupling_y;

	cavity_type cavity = cavity_from_name(p.get<std::string>("cavity"));
	double g = p.get<double>("g");
	if(cavity == cavity_type::polarized) {
		double phi = p.get<double>("phi");
		mode_coupling_x = {g * cos(phi)};
		mode_coupling_y = {g * sin(phi)};
	} else if(cavity == cavity_type::degenerate) {
		mode_coupling_x = {g, 0};
		mode_coupling_y = {0, g};
	}

	return {unitcells::square, {{spin_dim}}, {{J, mode_coupling_x}, {J, mode_coupling_y}}};
}

std::unique_ptr<cavity_magnet> cavity_magnet_from_param(const loadl::parser &p) {
	auto lat = p.get<std::string>("lattice");

	cavity_type cavity = cavity_from_name(p.get<std::string>("cavity"));
	std::vector<cavity_magnet::mode> modes;
	int max_photons = p.get<int>("max_photons");
	double omega = p.get<double>("omega");

	if(cavity == cavity_type::polarized) {
		modes = {{omega, max_photons}};
	} else if(cavity == cavity_type::degenerate) {
		modes = {{omega, max_photons}, {omega, max_photons}};
	}

	cavity_magnet_proto proto;
	int Lx = p.get<int>("Lx");
	int Ly = p.get<int>("Ly", Lx);

	double U = p.get<double>("U");

	if(lat == "square") {
		proto = make_square(p);
	}

	return std::make_unique<cavity_magnet>(lattice{proto.uc, Lx, Ly}, modes, proto.sites,
	                                       proto.bonds, U, cavity_magnet_measurement_settings{p});
}
