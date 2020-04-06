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
				throw std::runtime_error("lattice not bipartite!");
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

lattice::lattice(const unitcell &uc, int Lx, int Ly)
    : Lx{Lx}, Ly{Ly}  {
	int uc_spin_count = uc.sites.size();

	auto split_idx = [&](int i) { 
		return std::tuple{i%uc_spin_count, (i/uc_spin_count)%Lx, (i/uc_spin_count)/Lx};
	};

	for(int y = 0; y < Ly; y++) {
		for(int x = 0; x < Lx; x++) {
			for(auto &s : uc.sites) {
				sites.emplace_back(site{(x + s.pos[0]) * uc.a1 + (y + s.pos[1]) * uc.a2, s.nspinhalfs, s.Jin});
				spinhalf_count += s.nspinhalfs;
			}

			for(auto b : uc.bonds) {
				auto [iuc, ix, iy] = split_idx(b.i);
				auto [juc, jx, jy] = split_idx(b.j);

				int i =
				    Lx * uc_spin_count * ((iy + y) % Ly) + uc_spin_count * ((ix + x) % Lx) + iuc;
				int j =
				    Lx * uc_spin_count * ((jy + y) % Ly) + uc_spin_count * ((jx + x) % Lx) + juc;

				bonds.emplace_back(bond{i, j, b.J});
			}
		}
	}

	init_sublattice();
	init_vertex_data(uc);
}

void lattice::init_vertex_data(const unitcell &uc) {
	std::transform(uc.bonds.begin(), uc.bonds.end(), std::back_inserter(vertices_), [&](const bond &b) {
		return vertex_data{b, sites[b.i], sites[b.j]};
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
		vd.print(sites[b.i], sites[b.j]);
		idx++;
	}
}

void lattice::to_json(nlohmann::json &out) {
	for(const auto &site : sites) {
		out["sites"].push_back({
			{"pos", site.pos},
			{"sublattice", site.sublattice},
			{"nspinhalfs", site.nspinhalfs},
			{"Jin", site.Jin},
		});
	}

	for(const auto &bond : bonds) {
		out["bonds"].push_back({
		    {"i", bond.i},
		    {"j", bond.j},
		    {"J", bond.J},
		});
	}
}
