#include "cavity_magnet.h"
#include "nlohmann/json.hpp"
#include "util/kronecker_product.h"

cavity_magnet::cavity_magnet(const lattice &lat, const std::vector<mode>& modes, const std::vector<bond>& bonds) : model{model_type::cavity_magnet}, lat{lat}, modes_{modes}, bonds_{bonds} {
	assert(bonds_.size() == lat.uc.bonds.size());
}

// clang-format off
static const Eigen::MatrixXd heisenberg_exchange = init_mat(4,4, {
	1, 0, 0, 0,
	0,-1, 2, 0,
	0, 2,-1, 0,
	0, 0, 0, 1,
})/4;

static const Eigen::MatrixXd Szi = init_mat(4,4, {
	1, 0, 0, 0,
	0, 1, 0, 0,
	0, 0,-1, 0,
	0, 0, 0,-1,
})/2;

static const Eigen::MatrixXd Szj = init_mat(4,4, {
	1, 0, 0, 0,
	0,-1, 0, 0,
	0, 0, 1, 0,
	0, 0, 0,-1,
})/2;
// clang-format on
						
static Eigen::MatrixXd bond_hamiltonian(const cavity_magnet::mode& mode, const cavity_magnet::bond& bond, const unitcell::site& si, const unitcell::site& sj) {
	int dim = 4 * mode.max_bosons;
	Eigen::MatrixXd res = Eigen::MatrixXd::Zero(dim, dim);

	Eigen::MatrixXd boson_number = res.Zero(mode.max_bosons, mode.max_bosons);
	int n = 0;
	for(auto &d : boson_number.diagonal()) {
		d = n;
		n++;
	}

	return kronecker_prod(boson_number, Szi/si.coordination+Szj/sj.coordination) +
	bond.J * kronecker_prod(Eigen::MatrixXd::Identity(mode.max_bosons, mode.max_bosons), heisenberg_exchange) + mode.omega * kronecker_prod(boson_number, Eigen::MatrixXd::Identity(4,4));
	
}

sse_data cavity_magnet::generate_sse_data() const {
	std::vector<vertex_data> vert_data;
	std::vector<sse_data::bond> sse_bonds;
	std::vector<sse_data::site> sse_sites;

	int mode_count = modes_.size();

	int bond_idx{};
	for(const auto& b : bonds_) {
		for(const auto &m : modes_) {
			auto H = bond_hamiltonian(m, b, lat.uc.sites[lat.uc.bonds[bond_idx].i], lat.uc.sites[lat.uc.bonds[bond_idx].j.uc]);
			vert_data.push_back({{m.max_bosons, 2, 2}, H});
		}
		bond_idx++;
	}

	bond_idx = 0;
	for(const auto& b : lat.bonds) {
		for(int mode_idx = 0; mode_idx < static_cast<int>(modes_.size()); mode_idx++) {
			int bond_type = b.type*mode_count + mode_idx;
			sse_bonds.push_back({bond_type, {mode_idx, mode_count + b.i, mode_count + b.j}});
		}
		bond_idx++;
	}

	std::transform(modes_.begin(), modes_.end(), std::back_inserter(sse_sites), [](const auto& m) {
		return sse_data::site{m.max_bosons};
	});

	sse_sites.insert(sse_sites.end(), lat.Lx*lat.Ly*lat.uc.sites.size(), {2});

	return sse_data{vert_data, sse_sites, sse_bonds};
}

void cavity_magnet::to_json(nlohmann::json& out) const {
	out["model"] = "cavity_magnet";
	lat.to_json(out);
	for(const auto& m : modes_) {
		out["modes"].push_back({{"omega", m.omega}, {"coupling", m.coupling}, {"max_bosons", m.max_bosons}});
	}
}
