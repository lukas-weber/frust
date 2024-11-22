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

static std::vector<std::vector<double>> make_cavity(const unitcell &uc, const loadl::parser &p) {
	cavity_type cavity = cavity_from_name(p.get<std::string>("cavity"));
	double g = p.get<double>("g");

	std::vector<vec2> polarizations;
	if(cavity == cavity_type::polarized) {
		double phi = p.get<double>("phi");
		polarizations = {{cos(phi), sin(phi)}};
	} else if(cavity == cavity_type::degenerate) {
		polarizations = {{1, 0}, {0, 1}};
	}

	std::vector<std::vector<double>> mode_couplings;

	Eigen::Matrix2d A;
	A.col(0) = uc.a1;
	A.col(1) = uc.a2;
	for(const auto &bond : uc.bonds) {
		vec2 d = A * (vec2{bond.j.dx, bond.j.dy} + uc.sites[bond.j.uc].pos - uc.sites[bond.i].pos);
		mode_couplings.push_back({});
		for(const auto &pol : polarizations) {
			mode_couplings.back().push_back(g * pol.dot(d));
		}
	}
	return mode_couplings;
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

	const unitcell uc = unitcells::square;
	auto mode_couplings = make_cavity(uc, p);

	return {unitcells::square, {{spin_dim}}, {{J, mode_couplings[0]}, {J, mode_couplings[1]}}};
}

static cavity_magnet_proto make_honeycomb(const loadl::parser &p) {
	double J = p.get<double>("J");
	int spin_dim = round(2 * p.get<double>("S", 0.5) + 1);

	std::vector<double> mode_coupling_x;
	std::vector<double> mode_coupling_y;
	const auto &uc = unitcells::honeycomb;

	std::vector<std::vector<double>> mode_couplings = make_cavity(uc, p);

	std::vector<cavity_magnet::bond> bonds(uc.bonds.size());
	std::transform(mode_couplings.begin(), mode_couplings.end(), bonds.begin(),
	               [J](const auto &mc) {
		               return cavity_magnet::bond{J, mc};
	               });

	return {unitcells::honeycomb, {{spin_dim}, {spin_dim}}, bonds};
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
	} else if(lat == "honeycomb") {
		proto = make_honeycomb(p);
	} else {
		throw std::runtime_error{fmt::format("unknown lattice '{}'", lat)};
	}

	return std::make_unique<cavity_magnet>(lattice{proto.uc, Lx, Ly}, modes, proto.sites,
	                                       proto.bonds, U, cavity_magnet_measurement_settings{p});
}
