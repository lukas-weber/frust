#pragma once

#include "vertex_data.h"
#include <vector>
#include <iostream>

class sse_data {
public:
	struct site {
		int dim{};
	};
	using bond = std::vector<int>;
private:
	int bond_type_count_{};
	int site_type_count_{};
	std::vector<vertex_data> vertex_data_;
	std::vector<site> site_data_;
	
	std::vector<int> bonds_; // [bond_count; leg_count/2]
public:
	int nlegs{};
	int site_count{};
	double energy_offset{};

	sse_data(const std::vector<vertex_data> &vert_data, const std::vector<site>& site_data, const std::vector<bond> bonds)
		: bond_type_count_{static_cast<int>(vert_data.size())}, site_type_count_{static_cast<int>(site_data.size())}, vertex_data_{vert_data}, site_data_{site_data} {

		for(const auto &b : bonds) {
			if(static_cast<int>(b.size()) > nlegs/2) {
				nlegs = 2 * b.size();
			}
		}

		bonds_.resize(bonds.size()*nlegs/2, -1);

		int i{};
		for(const auto &b : bonds) {
			int j{};
			for(int s : b) {
				if(site_count <= s) {
					site_count = s + 1;
				}
				
				bonds_[i*nlegs/2 + j] = s;
				j++;
			}
			i++;
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

	const int *get_bond(int bond_idx) const { // [nlegs/2] // XXX: replace by std::span or something
		return &bonds_[(nlegs/2) * bond_idx]; 
	}

	int bond_count() const {
		return bonds_.size()/(nlegs/2);
	}

	void print() const {
		//int idx{};
		for(const auto& vd : vertex_data_) {
			//std::cout << fmt::format("\nbond {}-{}\n", bonds[idx].i, bonds[idx].j);
			vd.print();
		}
	}
};
