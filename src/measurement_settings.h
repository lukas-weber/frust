#pragma once

#include <fmt/format.h>
#include <loadleveller/parser.h>

struct measurement_settings {
	measurement_settings(const loadl::parser &p) {
		for(const auto &obs : p.get<std::vector<std::string>>("measure")) {
			if(obs == "j") {
				measure_j = true;
			} else if(obs == "jcorrlen") {
				measure_jcorrlen = true;
			} else if(obs == "mag") {
				measure_mag = true;
			} else if(obs == "sxmag") {
				measure_sxmag = true;
			} else if(obs == "symag") {
				measure_symag = true;
			} else if(obs == "sxsymag") {
				measure_sxsymag = true;
			} else if(obs == "sxsucmag") {
				measure_sxsucmag = true;
			} else if(obs == "chirality") {
				measure_chirality = true;
			} else {
				throw std::runtime_error{fmt::format("unknown observable '{}'", obs)};
			}
		}

		if(p.get<bool>("loopcorr_as_strucfac", false)) {
			loopcorr_as_strucfac = true;
		}

		if(!measure_j && measure_jcorrlen) {
			throw std::runtime_error{"jcorrlen can only be measured if j is measured"};
		}
		
		if(measure_sxsucmag && !p.defined("uc_signs")) {
			throw std::runtime_error{"wanted to measure unitcell staggered magnetization, but uc_signs not defined"};
		}
	}

	bool measure_j{};
	bool measure_mag{};
	bool measure_sxmag{};
	bool measure_symag{};
	bool measure_sxsymag{};
	bool measure_sxsucmag{};
	bool measure_chirality{};
	bool measure_jcorrlen{};

	bool loopcorr_as_strucfac{};
};
