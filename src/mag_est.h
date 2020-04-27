#pragma once

#include <loadleveller/measurements.h>
#include <loadleveller/evalable.h>
#include "lattice.h"
#include "opercode.h"

template <bool Staggered>
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
	double sign_{};
public:
	mag_est(const lattice &lat, double T, double sign) : lat_{lat}, T_{T}, sign_{sign} {
	}

	void init(const std::vector<state_idx> &spin) {
		tmpmag_ = 0;

		int i = 0;
		for(const auto &site : lat_.sites) {
			int sign = 1;
			if(Staggered) {
				sign *= site.sublattice;
			}
			tmpmag_ += sign * lat_.get_uc_site(i).basis.states[spin[i]].m;
			i++;
		}
		
		mag_ = tmpmag_;
		absmag_ = fabs(mag_);
		mag2_ = mag_ * mag_;
		mag4_ = mag2_ * mag2_;
	}

	void measure(opercode op, const std::vector<state_idx> &) {
		if(!Staggered) {
			return;
		}
		
		if(!op.diagonal()) {
			const auto &bond = lat_.bonds[op.bond()];
			const auto &bi = lat_.get_uc_site(bond.i).basis;
			const auto &bj = lat_.get_uc_site(bond.j).basis;

			double m20 = bi.states[op.leg_state(2)].m-bi.states[op.leg_state(0)].m;
			double m31 = bj.states[op.leg_state(3)].m - bj.states[op.leg_state(1)].m;
			
			tmpmag_ += lat_.sites[bond.i].sublattice*m20 + lat_.sites[bond.j].sublattice*m31;
		}

		mag_ += tmpmag_;
		absmag_ += fabs(tmpmag_);
		mag2_ += tmpmag_ * tmpmag_;
		mag4_ += tmpmag_ * tmpmag_ * tmpmag_ * tmpmag_;
		n_++;
	}

	void result(loadl::measurements &measure) {
		std::string p = "Sign";
		if(Staggered) {
			p += "Stag";
		}

		double norm = 1. / lat_.spinhalf_count;
		mag_ *= norm;
		absmag_ *= norm;
		mag2_ *= norm * norm;
		mag4_ *= norm * norm * norm * norm;

		measure.add(p + "Mag", sign_*mag_ / n_);
		measure.add(p + "AbsMag", sign_*absmag_ / n_);
		measure.add(p + "Mag2", sign_*mag2_ / n_);
		measure.add(p + "Mag4", sign_*mag4_ / n_);

		double chi = 1 / T_ / (n_ + 1) / n_ * (mag_ * mag_ + mag2_) * lat_.spinhalf_count;
		measure.add(p + "MagChi", sign_*chi);
	}

	void register_evalables(loadl::evaluator &eval) {
		std::string p = "Sign";
		std::string p2 = Staggered ? "Stag" : "";
		auto unsign = [](const std::vector<std::vector<double>> &obs) {
			return std::vector<double>{obs[0][0]/obs[1][0]};
		};
		
		eval.evaluate(p2 + "Mag", {p+p2+"Mag", "Sign"}, unsign);
		eval.evaluate(p2 + "Mag2", {p+p2+"Mag2", "Sign"}, unsign);
		eval.evaluate(p2 + "Mag4", {p+p2+"Mag4", "Sign"}, unsign);
		
		eval.evaluate(p2 + "MagChi", {p+p2+"MagChi", "Sign"}, unsign);
	}
};
