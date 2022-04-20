#pragma once

#include "../common/mag_est.h"
#include "cluster_magnet.h"
#include "j_est.h"
#include <tuple>

inline auto cluster_magnet_opstring_estimators(const cluster_magnet &cm, double T, double sign) {
	using M = cluster_magnet;
	auto obs = std::tuple{
	    j_est{cm, sign, cm.settings.measure_jcorrlen},
	    mag_est<mag_sign::none, M>{cm, T, sign},
	    mag_est<mag_sign::x, M>{cm, T, sign},
	    mag_est<mag_sign::y, M>{cm, T, sign},
	    mag_est<mag_sign::x | mag_sign::y, M>{cm, T, sign},
	    mag_est<mag_sign::x | mag_sign::uc, M>{cm, T, sign},
	};

	std::array flags = {
	    cm.settings.measure_j || cm.settings.measure_chirality,
	    cm.settings.measure_mag,
	    cm.settings.measure_sxmag,
	    cm.settings.measure_symag,
	    cm.settings.measure_sxsymag,
	    cm.settings.measure_sxsucmag,
	};
	return std::make_tuple(obs, flags);
}
