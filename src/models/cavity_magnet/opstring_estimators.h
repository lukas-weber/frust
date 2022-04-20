#pragma once

#include "../common/mag_est.h"
#include "cavity_magnet.h"
#include "photon_est.h"

auto cavity_magnet_opstring_estimators(const cavity_magnet &cm, double T, double sign) {
	using M = cavity_magnet;

	auto obs = std::tuple{
	    photon_est{cm, sign},
	    mag_est<mag_sign::none, M>{cm, T, sign},
	    mag_est<mag_sign::x | mag_sign::y, M>{cm, T, sign},
	    mag_est<mag_sign::x | mag_sign::uc, M>{cm, T, sign},
	};

	std::array flags = {
	    true,
	    cm.settings.measure_mag,
	    cm.settings.measure_sxsymag,
	    cm.settings.measure_sxsucmag,
	};

	return std::make_tuple(obs, flags);
}
