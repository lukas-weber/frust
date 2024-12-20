#pragma once

#include <loadleveller/parser.h>

struct cavity_magnet_measurement_settings {
	cavity_magnet_measurement_settings(const loadl::parser &p) {
		for(const auto &obs : p.get<std::vector<std::string>>("measure")) {
			if(obs == "mag") {
				measure_mag = true;
			} else if(obs == "sxsymag") {
				measure_sxsymag = true;
			} else if(obs == "sxsucmag") {
				measure_sxsucmag = true;
			} else if(obs == "sucmag") {
				measure_sucmag = true;
			} else {
				throw std::runtime_error{fmt::format("unknown observable '{}'", obs)};
			}
		}
	}

	bool measure_mag{};
	bool measure_sucmag{};
	bool measure_sxsymag{};
	bool measure_sxsucmag{};
};
