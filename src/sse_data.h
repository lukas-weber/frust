#pragma once

#include "vertex_data.h"
#include <vector>

class sse_data {
public:
	struct site {
		int dim{};
	};
	struct bond {
		int i{};
		int j{};
	};

private:
	int bond_type_count_{};
	int site_type_count_{};
	std::vector<vertex_data> vertex_data_;
	std::vector<site> site_data_;
public:
	int site_count{};
	std::vector<bond> bonds;
	double energy_offset{};

	sse_data(const std::vector<vertex_data> &vert_data, const std::vector<site>& site_data, const std::vector<bond> &bonds)
		: bond_type_count_{static_cast<int>(vert_data.size())}, site_type_count_{static_cast<int>(site_data.size())}, vertex_data_{vert_data}, site_data_{site_data}, bonds(bonds) {
			
		for(const auto &b : bonds) {
			if(site_count <= b.i) {
				site_count = b.i + 1;
			}
			if(site_count <= b.j) {
				site_count = b.j + 1;
			}
		}

		assert(site_count % site_type_count_ == 0);
		assert(bonds.size() % bond_type_count_ == 0);
		for(const auto &vd : vert_data) {
			energy_offset += vd.energy_offset;
		}
		energy_offset *= bonds.size()/bond_type_count_;
	}

	const vertex_data& get_vertex_data(int bond_idx) const {
		return vertex_data_[bond_idx % bond_type_count_];
	}

	const site& get_site_data(int site_idx) const {
		return site_data_[site_idx % site_type_count_];
	}
};
