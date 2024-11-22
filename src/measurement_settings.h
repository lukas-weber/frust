#pragma once

#include <loadleveller/parser.h>
#include <fmt/format.h>

struct measurement_settings {
	measurement_settings(const loadl::parser &p) {
		for(const auto &obs : p.get<std::vector<std::string>>("measure")) {
			if(obs == "j") {
				measure_j = true;
			} else if(obs == "mag") {
				measure_mag = true;
			} else if(obs == "sxmag") {
				measure_sxmag = true;
			} else if(obs == "symag") {
				measure_symag = true;
			} else if(obs == "sxsymag") {
				measure_sxsymag = true;
			} else {
				throw std::runtime_error{fmt::format("unknown observable '{}'", obs)};
			}
		}
	}
	
	bool measure_j{};
	bool measure_mag{};
	bool measure_sxmag{};
	bool measure_symag{};
	bool measure_sxsymag{};
};
