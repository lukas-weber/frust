#include "latticedef.h"

static unitcell make_square(int Lx, const loadl::parser &p) {
	unitcell uc;

	uc.a1 = {1,0};
	uc.a2 = {0,1};

	int nspinhalf{};
	std::string spin = p.get<std::string>("spin");
	if(spin == "1/2") {
		nspinhalf = 1;
	} else if(spin == "1") {
		nspinhalf = 2;
	} else {
		throw std::runtime_error{"illegal value for spin"};
	}

	double J = p.get<double>("J");

	uc.sites = {
		{{0,0}, nspinhalf, 0}
	};

	if(nspinhalf == 1) {
		uc.bonds = {
			{0,1,{J}},
			{0,Lx,{J}}
		};
	} else {
		uc.bonds = {
			{0,1,{J,J,J,J}},
			{0,Lx,{J,J,J,J}}
		};
	}

	return uc;
}

static unitcell make_triangle_square(int Lx, const loadl::parser &p) {
	unitcell uc;

	uc.a1 = {1,0};
	uc.a2 = {0,2};

	double Jtri = p.get<double>("Jtri");
	double Jin = p.get<double>("Jin",Jtri);

	double Jup = p.get<double>("Jup");

	double Jdirect = p.get<double>("Jdirect");
	double Jindirect = p.get<double>("Jindirect");

	uc.sites = {
		{{0,0}, 2, Jin},
		{{0,0.5}, 1, 0},
	};

	uc.bonds = {
		{0,1,{Jtri, Jtri}},
		{0,2,{Jindirect, Jindirect, Jdirect, Jindirect}},
		{1,2*Lx,{Jup, Jup}},
	};

	return uc;
}

lattice lattice_from_param(const loadl::parser &p) {
	auto name = p.get<std::string>("lattice");

	unitcell uc;

	int Lx = p.get<int>("Lx");
	int Ly = p.get<int>("Ly", Lx);

	if(name == "triangle_square") {
		uc = make_triangle_square(Lx, p);
	} else if(name == "square") {
		uc = make_square(Lx, p);
	} else {
		throw std::runtime_error{fmt::format("unknown lattice '{}'", name)};
	}

	return lattice{uc, Lx, Ly};
}
