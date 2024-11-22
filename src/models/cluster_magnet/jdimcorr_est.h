#pragma once

#include "lattice.h"
#include "opercode.h"
#include <loadleveller/evalable.h>
#include <loadleveller/measurements.h>

class jdimcorr_est {
private:
	const bool as_strucfac_{};
	const lattice &lat_;
	const double sign_{};

	double n_{1};
	double tmpjdim_{};

	std::vector<double> jdimcorr_{};

	int corrfunc_idx(int i, int j) const {
		if(as_strucfac_) {
			return 0;
		}
		auto [ix, iy, iuc] = lat_.split_idx(i);
		auto [jx, jy, juc] = lat_.split_idx(j);

		return lat_.Ly *
		           (lat_.Lx * (lat_.uc.sites.size() * iuc + juc) + (iy - jy + lat_.Ly) % lat_.Ly) +
		       (ix - jy + lat_.Lx) % lat_.Lx;
	}

public:
	jdimcorr_est(const lattice &lat, double sign, bool as_strucfac)
	    : as_strucfac_{as_strucfac}, lat_{lat}, sign_{sign} {
		if(!as_strucfac) {
			tmpjdimcorr_.resize(lat.uc.sites.size() * lat.sites.size(), 0);
			jdimcorr_.resize(tmpjdimcorr_.size(), 0);
		}
	}

	void init(const std::vector<state_idx> &spin) {
		for(int i = 0; i < lat_.sites.size(); i++) {
			for(int j = 0; j < lat_.sites.size(); j++) {
				int idx = corrfunc_idx(i, j);
				jdimcorr_[idx] += spin[i] * spin[j];
			}
		}
		tmpjdimcorr_ = jdimcorr_;
	}

	void measure(opercode op, const std::vector<state_idx> &spin) {
		if(!as_strucfac_) {
			return; // not implemented
		}
		int uc_size = lat_.uc.sites.size();
		if(!op.diagonal()) {
			const auto &bond = lat_.bonds[op.bond()];
			const auto &bi = lat_.get_uc_site(bond.i).basis;
			const auto &bj = lat_.get_uc_site(bond.j).basis;

			const auto &leg_state = lat_.get_vertex_data(op.bond()).get_legstate(op.vertex());

			double jdim20 = bi.states[leg_state[2]].jdim - bi.states[leg_state[0]].jdim;
			double jdim31 = bj.states[leg_state[3]].jdim - bj.states[leg_state[1]].jdim;

			tmp
		}
	}
};
