#pragma once

#include <loadleveller/measurements.h>
#include <loadleveller/evalable.h>
#include "lattice.h"
#include "opercode.h"

class mag_est {
private:
	double n_{1};

	double tmpmag_{};
	double mag_{};
	double absmag_{};
	double mag2_{};
	double mag4_{};

	const lattice &lat_;
	double T_{};
public:
	mag_est(const lattice &lat, double T) : lat_{lat}, T_{T} {
	}

	void init(const std::vector<jm> &spin) {
		tmpmag_ = 0;

		int i = 0;
		for(const auto &site : lat_.sites) {
			tmpmag_ += spin[i].m(site);
			i++;
		}
		
		mag_ = tmpmag_;
		absmag_ = fabs(mag_);
		mag2_ = mag_ * mag_;
		mag4_ = mag2_ * mag2_;
	}

	void measure(opercode op, const std::vector<jm> &) {
		return;
		if(!op.diagonal()) {
			const auto &bond = lat_.bonds[op.bond()];
			const auto &si = lat_.sites[bond.i];
			const auto &sj = lat_.sites[bond.j];

			tmpmag_ += (op.leg_state(2).m(si)-op.leg_state(0).m(si) + op.leg_state(3).m(sj) - op.leg_state(1).m(sj));
		}

		mag_ += tmpmag_;
		absmag_ += fabs(tmpmag_);
		mag2_ += tmpmag_ * tmpmag_;
		mag4_ += tmpmag_ * tmpmag_ * tmpmag_ * tmpmag_;
		n_++;
	}

	void result(loadl::measurements &measure) {
		std::string p = "";

		double norm = 1. / lat_.spinhalf_count;
		mag_ *= norm;
		absmag_ *= norm;
		mag2_ *= norm * norm;
		mag4_ *= norm * norm * norm * norm;

		measure.add(p + "Mag", mag_ / n_);
		measure.add(p + "AbsMag", absmag_ / n_);
		measure.add(p + "Mag2", mag2_ / n_);
		measure.add(p + "Mag4", mag4_ / n_);

		double chi = 1 / T_ / (n_ + 1) / n_ * (mag_ * mag_ + mag2_) * lat_.spinhalf_count;
		measure.add(p + "MagChi", chi);
	}

	void register_evalables(loadl::evaluator &) {

	}
};
