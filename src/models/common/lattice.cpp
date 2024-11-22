#include "lattice.h"
#include "nlohmann/json.hpp"

vec2 lattice::site_pos(int site_idx) const {
	auto [x, y, iuc] = split_idx(site_idx);
	return (x + uc.sites[iuc].pos[0]) * uc.a1 + (y + uc.sites[iuc].pos[1]) * uc.a2;
}

lattice::lattice(const unitcell &ucell, int Lx, int Ly) : uc{ucell}, Lx{Lx}, Ly{Ly} {
	int uc_spin_count = uc.sites.size();

	for(const auto &b : uc.bonds) {
		uc.sites[b.i].coordination++;
		uc.sites[b.j.uc].coordination++;
	}

	for(int y = 0; y < Ly; y++) {
		for(int x = 0; x < Lx; x++) {
			int bond_type{};
			for(auto b : uc.bonds) {
				assert(b.i < uc_spin_count);
				assert(b.j.uc < uc_spin_count);

				int i = Lx * uc_spin_count * (y % Ly) + uc_spin_count * (x % Lx) + b.i;
				int j = Lx * uc_spin_count * ((b.j.dy + y) % Ly) +
				        uc_spin_count * ((b.j.dx + x) % Lx) + b.j.uc;

				bonds.push_back(bond{bond_type, i, j});
				bond_type++;
			}
		}
	}
}

void lattice::to_json(nlohmann::json &out) const {
	out["Lx"] = Lx;
	out["Ly"] = Ly;
	out["uc_site_count"] = uc.sites.size();

	for(int idx = 0; idx < static_cast<int>(Lx * Ly * uc.sites.size()); idx++) {
		out["sites"].push_back({
		    {"pos", site_pos(idx)},
		});
	}
	for(const auto &bond : bonds) {
		out["bonds"].push_back({
		    {"i", bond.i},
		    {"j", bond.j},
		});
	}
}
