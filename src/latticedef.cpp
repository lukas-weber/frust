#include "latticedef.h"

static unitcell make_square(const loadl::parser &p) {
	unitcell uc;

	uc.a1 = {1,0};
	uc.a2 = {0,1};

	double J = p.get<double>("J");

	uc.sites = {
		{{0,0}, {}, site_bases::spin}
	};
	uc.bonds = {
		{0,{1,0,0},{J}},
		{0,{0,1,0},{J}}
	};

	return uc;
}

static unitcell make_bilayer(const loadl::parser &p) {
	unitcell uc;

	uc.a1 = {1,0};
	uc.a2 = {0,1};

	double Jpar = p.get<double>("Jpar");
	double Jperp = p.get<double>("Jperp");
	double Jx = p.get<double>("Jx", Jpar);

	uc.sites = {
		{{0,0}, {Jperp}, site_bases::dimer}
	};
	uc.bonds = {
		{0,{1,0,0},{Jpar,Jx,Jx,Jpar}},
		{0,{0,1,0},{Jpar,Jx,Jx,Jpar}}
	};

	return uc;
}

static unitcell make_dimerized_bilayer(const loadl::parser &p) {
	unitcell uc;

	uc.a1 = {1,0};
	uc.a2 = {0,2};

	double Jpar = p.get<double>("Jpar");
	double Jperp = p.get<double>("Jperp");
	double Jpardim = p.get<double>("Jpardim");

	uc.sites = {
		{{0,0}, {Jperp}, site_bases::dimer},
		{{0,0.5}, {Jperp}, site_bases::dimer}
	};
	uc.bonds = {
		{0,{0,0,1},{Jpardim,Jpardim,Jpardim,Jpardim}},
		{0,{1,0,0},{Jpar,Jpar,Jpar,Jpar}},
		{1,{1,0,1},{Jpar,Jpar,Jpar,Jpar}},
		{1,{0,1,0},{Jpar,Jpar,Jpar,Jpar}}
	};

	return uc;
}

static unitcell make_shastry_sutherland(const loadl::parser &p) {
	unitcell uc;
	double J = p.get<double>("J");
	double JD = p.get<double>("JD");

	uc.a1 = {1,0};
	uc.a2 = {0,1};

	std::string basis = p.get<std::string>("basis");

	if(basis == "spin") {
		uc.sites = {
			{{0,0}, {}, site_bases::spin},
			{{0.5,0}, {}, site_bases::spin},
			{{0.5,0.5}, {}, site_bases::spin},
			{{0,0.5}, {}, site_bases::spin},
		};

		uc.bonds = {
			{0,{0,0,1}, {J}},
			{1,{0,0,2}, {J}},
			{2,{0,0,3}, {J}},
			{3,{0,0,0}, {J}},
			{1,{0,0,3}, {JD}},
			
			{1,{1,0,0}, {J}},
			{2,{1,0,3}, {J}},
			{2,{0,1,1}, {J}},
			{3,{0,1,0}, {J}},
			{2,{1,1,0}, {JD}},
		};
	} else if(basis == "dimer") {
		uc.sites = {
			{{0,0}, {JD}, site_bases::dimer},
			{{0.5,0.5}, {JD}, site_bases::dimer},
		};

		uc.bonds = {
			{0, {0,0,1}, {0, J, 0, J}},
			{1, {1,0,0}, {0, J, 0, J}},
			{1, {0,1,0}, {J, 0, J, 0}},
			{1, {1,1,0}, {J, J, 0, 0}},
		};
	}



	//uc.sites = {
	//	{{-0.25,0}, {}, site_bases::spin},
	//	{{0,-0.25}, {}, site_bases::spin},
	//	{{0.25,0}, {}, site_bases::spin},
	//	{{0,0.25}, {}, site_bases::spin},
	//	{{-0.25+0.5,0.5}, {}, site_bases::spin},
	//	{{0.5,-0.25+0.5}, {}, site_bases::spin},
	//	{{0.25+0.5,0.5}, {}, site_bases::spin},
	//	{{0.5,0.25+0.5}, {}, site_bases::spin},
	//};

	//uc.bonds = {
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
	
	return uc;
}

static unitcell make_triangle(const loadl::parser &p) {
	unitcell uc;
	double J = p.get<double>("J");

	uc.a1 = {1,0};
	uc.a2 = {0,1};
	
	uc.sites = {
		{{0,0}, {}, site_bases::spin},
	};

	uc.bonds = {
		{0,{0,1,0}, {J}},
		{0,{1,0,0}, {J}},
		{0,{1,1,0}, {J}},
	};
	
	return uc;
}

static unitcell make_kagome(const loadl::parser &p) {
	unitcell uc;

	uc.a1 = {1,0};
	uc.a2 = {-0.5,sqrt(3)/2};
	
	double J3 = p.get<double>("J3");
	double J1 = p.get<double>("J1", J3);
	double J2 = p.get<double>("J2", J1);

	double J = p.get<double>("J");

	double h = p.get<double>("h",0.);

	std::string basis = p.get<std::string>("basis");

	if(basis == "trimer") {
		uc.sites = {
			{{0,0}, {J3,J1,J2}, site_bases::trimer,h},
		};

		uc.bonds = {
			{0,{1,0,0}, {0,0,0,
				     J,0,0,
				     0,0,0}},
			{0,{1,1,0}, {0,0,0,
				     0,0,0,
				     J,0,0,}},
			{0,{0,1,0}, {0,0,0,
				     0,0,0,
				     0,J,0}},
		};
	} else if(basis == "dimer") {
		uc.sites = {
			{{0,0}, {J3}, site_bases::dimer,h},
			{{0.2,0.5}, {}, site_bases::spin,h},
		};
		uc.bonds = {
			{0, {0,0,1}, {J1, J2}},
			{0, {1,0,0}, {0,J,
						  0,0}},
			{1, {0,1,0}, {J,0}},
			{1, {1,1,0}, {0,J}},
		};
	} else if(basis == "spin") {
		uc.sites = {
			{{0,0}, {}, site_bases::spin,h},
			{{0.5,0}, {}, site_bases::spin,h},
			{{0.5,0.5}, {}, site_bases::spin,h},
		};
		
		uc.bonds = {
			{0,{0,0,1}, {J3}},
			{1,{0,0,2}, {J1}},
			{2,{0,0,0}, {J2}},
			{1,{1,0,0}, {J}},
			{2,{1,1,0}, {J}},
			{2,{0,1,1}, {J}},
		};
	} else {
		throw std::runtime_error(fmt::format("unknown basis: {}", basis));
	}
	
	return uc;
}

static unitcell make_triangle_square(const loadl::parser &p) {
	unitcell uc;

	uc.a1 = {1,0};
	uc.a2 = {0,2};

	double J3 = p.get<double>("J3");
	double J1 = p.get<double>("J1", J3);
	double J2 = p.get<double>("J2", J1);

	double Jn = p.get<double>("Jn");
	double Jnn = p.get<double>("Jnn");

	int basis = p.get<int>("basis",2);

	if(basis == 1) {
		uc.sites = {
			{{-0.3,0}, {}, site_bases::spin},
			{{0.3,0}, {}, site_bases::spin},
			{{0,0.5}, {}, site_bases::spin},
		};
		uc.bonds = {
			{0,{0,0,1}, {J3}},
			{1,{0,0,2}, {J2}},
			{2,{0,0,0}, {J1}},
			{1,{1,0,0}, {Jn}},
			{2,{0,1,0}, {Jn}},
			{2,{0,1,1}, {Jn}},
		};
		if(Jnn != 0) {
			throw std::runtime_error{"basis = 1 does not support Jnn"};
		}
	} else if(basis == 2) {
		uc.sites = {
			{{0,0}, {J3}, site_bases::dimer},
			{{0,0.5}, {}, site_bases::spin},
		};

		uc.bonds = {
			{0,{0,0,1},{J1, J2}},
			{0,{1,0,0},{Jnn, Jn, Jnn, Jnn}},
			{1,{0,1,0},{Jn, Jn}},
		};
	} else if(basis == 5) {
		uc.sites = {
			{{0.5,0.5}, {J2}, site_bases::dimer},
			{{0,0}, {}, site_bases::spin},
		};

		uc.bonds = {
			{0,{0,0,1},{J3, J1}},
			{0,{1,0,1},{0, Jn}},
			{0,{0,1,0},{0, Jn, 0, 0}},
			{0,{0,1,1},{Jn, 0}},
		};
	} else if(basis == 3 || basis == 4) {
		uc.sites = {{{0,0}, {J3, J1, J2}, basis == 3 ? site_bases::trimer : site_bases::trimer23}};

		uc.bonds = {
			{0, {1,0,0}, {Jnn,Jnn,Jnn,
				      Jn, Jnn, Jnn,
				      Jnn, Jnn, Jnn}},
			{0, {0,1,0}, {Jnn, Jnn, Jnn,
				      Jnn, Jnn, Jnn,
				      Jn, Jn, Jnn}},
		};
	} else {
		throw std::runtime_error("invalid basis! must be 1, 2 or 3");
	}

	return uc;
}

lattice lattice_from_param(const loadl::parser &p, bool with_vertex_data) {
	auto name = p.get<std::string>("lattice");

	unitcell uc;

	int Lx = p.get<int>("Lx");
	int Ly = p.get<int>("Ly", Lx);

	if(name == "triangle_square") {
		uc = make_triangle_square(p);
	} else if(name == "square") {
		uc = make_square(p);
	} else if(name == "bilayer") {
		uc = make_bilayer(p);
	} else if(name == "dimerized_bilayer") {
		uc = make_dimerized_bilayer(p);
	} else if(name == "shastry_sutherland") {
		uc = make_shastry_sutherland(p);
	} else if(name == "triangle") {
		uc = make_triangle(p);
	} else if(name == "kagome") {
		uc = make_kagome(p);
	} else {
		throw std::runtime_error{fmt::format("unknown lattice '{}'", name)};
	}

	return lattice{uc, Lx, Ly, with_vertex_data};
}
