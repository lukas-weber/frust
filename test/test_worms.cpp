#include "worms.h"
#include <catch2/catch.hpp>

TEST_CASE("worms") {
	int dim = GENERATE(2, 4, 7, 8);
	std::vector<int> s_ins;
	for(state_idx s_in = 0; s_in < dim; s_in++) {
		s_ins.push_back(s_in);
	}
	
	for(int s_in : s_ins) {
		std::vector<int> s_outs{s_in};
		
		for(worm_idx w = 0; w < worm_count(dim); w++) {
			s_outs.push_back(worm_action(w, s_in, dim)); 
			REQUIRE(worm_action(worm_inverse(w, dim), s_outs.back(), dim) == s_in); 
		}
		std::sort(s_outs.begin(), s_outs.end());
		REQUIRE(s_ins == s_outs);
	}
}
		
		
