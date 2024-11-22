#include "lattice.h"
#include <tuple>

void lattice::init_sublattice() {
	if(sites.size() == 0) {
		return;
	}

	sublattice.resize(sites.size(), 0);
	sublattice.at(0) = 1;

	int num_set = 1;

	for(int tries = 0; tries < sites.size(); tries++) {
		for(const auto &b : bonds) {
			if(sublattice[b.i] != 0 && sublattice[b.j] == sublattice[b.i]) {
				throw std::runtime_error("lattice not bipartite!");
			}

			if(sublattice[b.i] != 0) {
				sublattice[b.j] = -sublattice[b.i];
				num_set++;
			}

			if(sublattice[b.j] != 0) {
				sublattice[b.i] = -sublattice[b.j];
				num_set++;
			}
		}
		if(num_set >= sites.size())
			break;
	}
}

lattice::lattice(unitcell uc, int Lx, int Ly)
    : Lx{Lx}, Ly{Ly}  {
	int uc_spin_count = uc.sites.size();

	auto split_idx = [&](int i) { 
		return std::tuple{i%uc_spin_count, (i/uc_spin_count)%Lx, (i/uc_spin_count)/Lx};
	};

	for(int y = 0; y < Ly; y++) {
		for(int x = 0; x < Lx; x++) {
			for(auto &s : uc.sites) {
				sites.emplace_back(site{(x + s.pos[0]) * uc.a1 + (y + s.pos[1]) * uc.a2, s.Jin});
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
}
