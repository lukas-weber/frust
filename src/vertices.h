#pragma once

#include "bond.h"
#include <vector>
#include "jm.h"

struct transition {
	std::array<double,4*4> probs;
	std::array<int,4*4> targets;
};

class vertex_data {
public: 
	static const int basis_size = 4;
	static const int leg_count = 4;
	
	std::vector<transition> transitions; // [vertex*4*basis_size + worm_action_in*4 + leg_in]
	
	vertex_data(const bond &b, const site &si, const site &sj);
	void print() const;
private:
	std::vector<std::array<jm,4>> legstates_;
	std::vector<double> weights_;
	
	void construct_vertices(const bond &b, const site &si, const site &sj);
	int vertex_change_apply(int vertex, int leg_in, jm_action action_in, int leg_out, jm_action action_out) const;
};
