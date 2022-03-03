#include "cluster_magnet_def.h"
#include "cluster_magnet.h"
#include "basis.h"

struct cluster_magnet_proto {
	unitcell uc;
	std::vector<cluster_site> sites;
	std::vector<cluster_bond> bonds;
};

static cluster_magnet_proto make_square(const loadl::parser &p) {
	double J = p.get<double>("J");

	return {unitcells::square,
		{{{}, site_bases::spin}},
		{{{J}},
		 {{J}}}};
}

static cluster_magnet_proto make_bilayer(const loadl::parser &p) {
	double Jpar = p.get<double>("Jpar");
	double Jperp = p.get<double>("Jperp");
	double Jx = p.get<double>("Jx", Jpar);

	return {unitcells::square,
		{{{Jperp}, site_bases::dimer}},
		{{{Jpar, Jx, Jx, Jpar}}, {{Jpar, Jx, Jx, Jpar}}}};
}

static cluster_magnet_proto make_dimerized_bilayer(const loadl::parser &p) {
	double Jpar = p.get<double>("Jpar");
	double Jperp = p.get<double>("Jperp");
	double Jpardim = p.get<double>("Jpardim");

	return {unitcells::columnar_dimer,
		{{{Jperp}, site_bases::dimer, 0, 1},
		 {{Jperp}, site_bases::dimer, 0, -1}},
		{{{Jpardim, Jpardim, Jpardim, Jpardim}},
		 {{Jpar, Jpar, Jpar, Jpar}},
		 {{Jpar, Jpar, Jpar, Jpar}},
		 {{Jpar, Jpar, Jpar, Jpar}}}};
}

static cluster_magnet_proto make_shastry_sutherland(const loadl::parser &p) {
	double J = p.get<double>("J");
	double JD = p.get<double>("JD");
	std::string basis = p.get<std::string>("basis");

	if(basis == "spin") {
		return {unitcells::shastry_sutherland,
			{{{}, site_bases::spin},
			 {{}, site_bases::spin},
			 {{}, site_bases::spin},
			 {{}, site_bases::spin}},
			{{{J}},
			 {{J}},
			 {{J}},
			 {{J}},
			 {{JD}},
			 {{J}},
			 {{J}},
			 {{J}},
			 {{J}},
			 {{JD}}}
		};
	}

	if(basis == "dimer") {
		return {unitcells::shastry_sutherland_dimer,
		       {{{JD}, site_bases::dimer},
		        {{JD}, site_bases::dimer}},
		       {{{0, J, 0, J}},
		        {{0, J, 0, J}},
		        {{J, 0, J, 0}},
		        {{J, J, 0, 0}}}
		};
	}

	// uc.sites = {
	//	{{-0.25,0}, {}, site_bases::spin},
	//	{{0,-0.25}, {}, site_bases::spin},
	//	{{0.25,0}, {}, site_bases::spin},
	//	{{0,0.25}, {}, site_bases::spin},
	//	{{-0.25+0.5,0.5}, {}, site_bases::spin},
	//	{{0.5,-0.25+0.5}, {}, site_bases::spin},
	//	{{0.25+0.5,0.5}, {}, site_bases::spin},
	//	{{0.5,0.25+0.5}, {}, site_bases::spin},
	//};

	// uc.bonds = {
	//	{0,{0,0,1}, {J}},
	//	{1,{0,0,2}, {J}},
	//	{2,{0,0,3}, {J}},
	//	{3,{0,0,0}, {J}},
	//	{1,{0,0,3}, {JD}},
	//
	//	{4,{0,0,5}, {J}},
	//	{5,{0,0,6}, {J}},
	//	{6,{0,0,7}, {J}},
	//	{7,{0,0,4}, {J}},
	//	{5,{0,0,7}, {JD}},
	//
	//	{2,{0,0,5}, {J}},
	//	{3,{0,0,4}, {J}},
	//
	//	{5,{1,0,0}, {J}},
	//	{6,{1,0,3}, {J}},
	//	{2,{1,0,0}, {JD}},
	//	{6,{1,0,4}, {JD}},
	//
	//	{4,{0,1,1}, {J}},
	//	{7,{0,1,2}, {J}},
	//	{6,{1,1,1}, {J}},
	//	{7,{1,1,0}, {J}},
	//};

	throw std::runtime_error{fmt::format("unknown basis {}", basis)};
}

static cluster_magnet_proto make_triangle(const loadl::parser &p) {
	double J = p.get<double>("J");
	return {unitcells::triangle, {{{}, site_bases::spin}}, {{{J}},{{J}},{{J}}}};
}

static cluster_magnet_proto make_kagome(const loadl::parser &p) {
	double J3 = p.get<double>("J3");
	double J1 = p.get<double>("J1", J3);
	double J2 = p.get<double>("J2", J1);

	double J = p.get<double>("J");

	double h = p.get<double>("h", 0.);

	std::string basis = p.get<std::string>("basis");

	if(basis == "trimer") {
		return {unitcells::kagome_trimer,
			{{{J3, J1, J2}, site_bases::trimer, h}},
		    	{{{0, 0, 0, 0, 0, 0, 0, J, 0}},
			 {{0, 0, 0, 0, 0, 0, J, 0, 0}},
			 {{0, 0, 0, J, 0, 0, 0, 0, 0}}}
		};
	}
	if(basis == "dimer") {
		return {unitcells::kagome_dimer_spin, 
			{{{J3}, site_bases::dimer, h},
		  	 {{}, site_bases::spin, h}},
		  	{{{J1, J2}},
		         {{0, J, 0, 0}},
		         {{J, 0}},
		         {{0, J}}}
		};
	}
	if(basis == "spin") {
		return {unitcells::kagome,
			{{{}, site_bases::spin, h},
		         {{}, site_bases::spin, h},
		         {{}, site_bases::spin, h}},
			{{{J3}}, {{J1}}, {{J2}},
			 {{J}},  {{J}},  {{J}}}
		};
	}
	
	throw std::runtime_error(fmt::format("unknown basis: {}", basis));
}

static cluster_magnet_proto make_diamond_square(const loadl::parser &p) {
	double J1 = p.get<double>("J1");
	double J2 = p.get<double>("J2", J1);

	double h = p.get<double>("h", 0.);

	return {unitcells::lieb_lattice,
		{{{}, site_bases::spin, h},
		 {{J2}, site_bases::dimer, h},
		 {{J2}, site_bases::dimer, h}},
		{{{J1, J1}},
		 {{J1, J1}},
		 {{J1, J1}},
		 {{J1, J1}}}
	};
}

static cluster_magnet_proto make_triangle_square(const loadl::parser &p) {
	double J3 = p.get<double>("J3");
	double J1 = p.get<double>("J1", J3);
	double J2 = p.get<double>("J2", J1);

	double Jn = p.get<double>("Jn");
	double Jnn = p.get<double>("Jnn");

	int basis = p.get<int>("basis", 2);

	if(basis == 1) {
		if(Jnn != 0) {
			throw std::runtime_error{"basis = 1 does not support Jnn"};
		}
		return {unitcells::triangle_square,
			{{{}, site_bases::spin},
			 {{}, site_bases::spin},
			 {{}, site_bases::spin}},
			{{{J3}}, {{J2}}, {{J1}},
			 {{Jn}}, {{Jn}}, {{Jn}}}
		};
	}
	if(basis == 2) {
		return {unitcells::triangle_square_dimer_spin,
		        {{{J3}, site_bases::dimer},
		         {{}, site_bases::spin}},
		        {{{J1, J2}},
		         {{Jnn, Jn, Jnn, Jnn}},
		         {{Jn, Jn}}}
		};
	}
	if(basis == 5) {
		return {unitcells::triangle_square_dimer_spin_rot,
			{{{J2}, site_bases::dimer},
			 {{}, site_bases::spin}},
		        {{{J3, J1}},
		         {{0, Jn}},
		         {{0, Jn, 0, 0}},
		         {{Jn, 0}}}
		};
	}
	if(basis == 3 || basis == 4) {
		return {unitcells::triangle_square_trimer,
			{{{J3, J1, J2}, basis == 3 ? site_bases::trimer : site_bases::trimer23}},
			{{{Jnn, Jnn, Jnn, Jn, Jnn, Jnn, Jnn, Jnn, Jnn}},
		         {{Jnn, Jnn, Jnn, Jnn, Jnn, Jnn, Jn, Jn, Jnn}}}
		};
	}

	throw std::runtime_error(fmt::format("unknown basis {}", basis));
}

std::unique_ptr<cluster_magnet> cluster_magnet_from_param(const loadl::parser &p) {
	auto name = p.get<std::string>("lattice");

	cluster_magnet_proto proto;

	int Lx = p.get<int>("Lx");
	int Ly = p.get<int>("Ly", Lx);

	if(name == "triangle_square") {
		proto = make_triangle_square(p);
	} else if(name == "square") {
		proto = make_square(p);
	} else if(name == "diamond_square") {
		proto = make_diamond_square(p);
	} else if(name == "bilayer") {
		proto = make_bilayer(p);
	} else if(name == "dimerized_bilayer") {
		proto = make_dimerized_bilayer(p);
	} else if(name == "shastry_sutherland") {
		proto = make_shastry_sutherland(p);
	} else if(name == "triangle") {
		proto = make_triangle(p);
	} else if(name == "kagome") {
		proto = make_kagome(p);
	} else {
		throw std::runtime_error{fmt::format("unknown lattice '{}'", name)};
	}

	return std::make_unique<cluster_magnet>(lattice{proto.uc, Lx, Ly}, proto.sites, proto.bonds);
}
