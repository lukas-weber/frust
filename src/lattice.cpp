#include "lattice.h"
#include <tuple>

lattice::lattice(const unitcell &ucell, int Lx, int Ly, bool with_vertex_data)
    : uc{ucell}, Lx{Lx}, Ly{Ly} {
	int uc_spin_count = uc.sites.size();

	if(Ly == 1) {
		auto it = std::remove_if(uc.bonds.begin(), uc.bonds.end(),
		                         [](const auto &b) { return b.j.dy != 0; });
		uc.bonds.erase(it, uc.bonds.end());
	}
	if(Lx == 1) {
		auto it = std::remove_if(uc.bonds.begin(), uc.bonds.end(),
		                         [](const auto &b) { return b.j.dx != 0; });
		uc.bonds.erase(it, uc.bonds.end());
	}

	for(const auto &b : uc.bonds) {
		uc.sites[b.i].coordination++;
		uc.sites[b.j.uc].coordination++;
	}

	for(int y = 0; y < Ly; y++) {
		for(int x = 0; x < Lx; x++) {
			for(auto &s : uc.sites) {
				sites.emplace_back(lat_site{(x + s.pos[0]) * uc.a1 + (y + s.pos[1]) * uc.a2});
				spinhalf_count += s.basis.nspinhalfs;
			}

			for(auto b : uc.bonds) {
				assert(b.i < uc_spin_count);
				assert(b.j.uc < uc_spin_count);

				int i = Lx * uc_spin_count * (y % Ly) + uc_spin_count * (x % Lx) + b.i;
				int j = Lx * uc_spin_count * ((b.j.dy + y) % Ly) +
				        uc_spin_count * ((b.j.dx + x) % Lx) + b.j.uc;

				bonds.emplace_back(lat_bond{i, j});
			}
		}
	}

	if(with_vertex_data) {
		init_vertex_data(uc);
	}
}

void lattice::init_vertex_data(const unitcell &uc) {
	std::transform(uc.bonds.begin(), uc.bonds.end(), std::back_inserter(vertices_),
	               [&](const uc_bond &b) {
		               return vertex_data{b, uc.sites[b.i], uc.sites[b.j.uc]};
	               });
	energy_offset = 0;
	for(const auto &vd : vertices_) {
		energy_offset += vd.energy_offset;
	}
	energy_offset *= Lx * Ly;
}

void lattice::vertex_print() const {
	int idx{};
	for(const auto &vd : vertices_) {
		const auto &b = bonds[idx];
		vd.print(get_uc_site(b.i).basis, get_uc_site(b.j).basis);
		idx++;
	}
}

void lattice::to_json(nlohmann::json &out) {
	int idx{};

	out["Lx"] = Lx;
	out["Ly"] = Ly;
	out["uc_spin_count"] = uc.sites.size();

	for(const auto &site : sites) {
		const auto &uc_st = uc.sites[idx % uc.sites.size()];
		out["sites"].push_back({
		    {"pos", site.pos},
		    {"nspinhalfs", uc_st.basis.nspinhalfs},
		    {"Jin", uc_st.Jin},
		    {"h", uc_st.h},
		});
		idx++;
	}

	idx = 0;
	for(const auto &bond : bonds) {
		out["bonds"].push_back({
		    {"i", bond.i},
		    {"j", bond.j},
		    {"J", uc.bonds[idx % uc.bonds.size()].J},
		});
		idx++;
	}
}
