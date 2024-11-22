#include "lattice.h"
#include <tuple>
	
void lattice::init_sublattice() {
	if(sites.size() == 0) {
		return;
	}

	sites.at(0).sublattice = 1;

	int num_set = 1;

	for(const auto &site : sites) {
		(void)site;
		for(const auto &b : bonds) {
			if(sites[b.i].sublattice != 0 && sites[b.j].sublattice == sites[b.i].sublattice) {
				return;
			}

			if(sites[b.i].sublattice != 0) {
				sites[b.j].sublattice = -sites[b.i].sublattice;
				num_set++;
			}

			if(sites[b.j].sublattice != 0) {
				sites[b.i].sublattice = -sites[b.j].sublattice;
				num_set++;
			}
		}
		if(num_set >= static_cast<int>(sites.size())) {
			break;
		}
	}
}

lattice::lattice(const unitcell &ucell, int Lx, int Ly)
    : uc{ucell}, Lx{Lx}, Ly{Ly}  {
	int uc_spin_count = uc.sites.size();

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

				int i =
				    Lx * uc_spin_count * (y % Ly) + uc_spin_count * (x % Lx) + b.i;
				int j =
				    Lx * uc_spin_count * ((b.j.dy + y) % Ly) + uc_spin_count * ((b.j.dx + x) % Lx) + b.j.uc;

				bonds.emplace_back(lat_bond{i, j});
			}
		}
	}

	init_sublattice();
	init_vertex_data(uc);
}

void lattice::init_vertex_data(const unitcell &uc) {
	std::transform(uc.bonds.begin(), uc.bonds.end(), std::back_inserter(vertices_), [&](const uc_bond &b) {
		return vertex_data{b, uc.sites[b.i], uc.sites[b.j.uc]};
	});
	energy_offset = 0;
	for(const auto &vd : vertices_) {
		energy_offset += vd.energy_offset;
	}
	energy_offset *= Lx*Ly;
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
	for(const auto &site : sites) {
		const auto &uc_st = uc.sites[idx%uc.sites.size()];
		out["sites"].push_back({
			{"pos", site.pos},
			{"sublattice", site.sublattice},
			{"nspinhalfs", uc_st.basis.nspinhalfs},
			{"Jin", uc_st.Jin},
		});
		idx++;
	}
	
	idx = 0;
	for(const auto &bond : bonds) {
		out["bonds"].push_back({
		    {"i", bond.i},
		    {"j", bond.j},
		    {"J", uc.bonds[idx%uc.bonds.size()].J},
		});
		idx++;
	}
}
