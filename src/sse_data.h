#pragma once

#include "vertex_data.h"
#include <vector>

class sse_data {
public:
	struct site {
		int dim{};
	};
	struct bond {
		int type{};
		std::vector<int> sites;
	};
	int nlegs{};

private:
	std::vector<vertex_data> vertex_data_;

	std::vector<site> sites_;
	std::vector<int> bonds_; // [bond_count; bonds_stride()]
	int bonds_stride() const {
		return nlegs / 2 + 1;
	}

public:
	int site_count{};
	int bond_count{};
	double energy_offset{};

	sse_data(const std::vector<vertex_data> &vert_data, const std::vector<site> &sites,
	         const std::vector<bond> &bonds)
	    : vertex_data_{vert_data}, sites_{sites}, site_count{static_cast<int>(sites.size())},
	      bond_count{static_cast<int>(bonds.size())} {
		for(const auto &vd : vertex_data_) {
			if(vd.leg_count > nlegs) {
				nlegs = vd.leg_count;
			}
		}

		bonds_.resize(bonds.size() * bonds_stride(), -1);

		int i{};
		for(const auto &b : bonds) {
			assert(static_cast<int>(b.sites.size()) == vertex_data_[b.type].leg_count / 2);

			int j{};
			bonds_[i * bonds_stride()] = b.type;
			for(int s : b.sites) {
				assert(vertex_data_[b.type].dims[j] == sites[s].dim);
				assert(s < site_count);

				bonds_[i * bonds_stride() + j + 1] = s;
				j++;
			}
			i++;
			energy_offset += vertex_data_[b.type].energy_offset;
		}
	}

	const vertex_data &get_vertex_data(int bond_idx) const {
		return vertex_data_[bonds_[bond_idx * bonds_stride()]];
	}

	const site &get_site_data(int site_idx) const {
		return sites_[site_idx];
	}

	const int *get_bond(int bond_idx) const { // [nlegs/2] // XXX: replace by std::span or something
		return &bonds_[bonds_stride() * bond_idx + 1];
	}

	void print() const {
		// int idx{};
		for(const auto &vd : vertex_data_) {
			// std::cout << fmt::format("\nbond {}-{}\n", bonds[idx].i, bonds[idx].j);
			vd.print();
		}
	}
};
