#include "lattice.h"

vec2 lattice::site_pos(int site_idx) const {
	auto [x, y, iuc] = split_idx(site_idx);
	return (x + uc.sites[iuc].pos[0]) * uc.a1 + (y + uc.sites[iuc].pos[1]) * uc.a2;
}
	

lattice::lattice(const unitcell &ucell, int Lx, int Ly)
    : uc{ucell}, Lx{Lx}, Ly{Ly}  {
	int uc_spin_count = uc.sites.size();
	// FIXME: replace by something more sensible
	/*if(Ly == 1) {
		auto it = std::remove_if(uc.bonds.begin(), uc.bonds.end(),
		                         [](const auto &b) { return b.j.dy != 0; });
		uc.bonds.erase(it, uc.bonds.end());
	}
	if(Lx == 1) {
		auto it = std::remove_if(uc.bonds.begin(), uc.bonds.end(),
		                         [](const auto &b) { return b.j.dx != 0; });
		uc.bonds.erase(it, uc.bonds.end());
	}*/

	for(const auto &b : uc.bonds) {
		uc.sites[b.i].coordination++;
		uc.sites[b.j.uc].coordination++;
	}

	for(int y = 0; y < Ly; y++) {
		for(int x = 0; x < Lx; x++) {
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
}
