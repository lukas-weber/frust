#pragma once

#include <loadleveller/measurements.h>
#include <loadleveller/evalable.h>
#include "lattice.h"
#include "opercode.h"

class j_est {
private:
	double n_{1};

	double tmpj_{};
	double j_{};
	double j2_{};


	const lattice &lat_;
	double sign_{};
public:
	j_est(const lattice &lat, double sign) : lat_{lat}, sign_{sign} {
	}

	void init(const std::vector<state_idx> &spin) {
		tmpj_ = 0;

		int i = 0;
		for(auto s : spin) {
			tmpj_ += lat_.get_uc_site(i).basis.states[s].j;
			i++;
		}
		
		j_ = tmpj_;
		j2_ = j_ * j_;
	}

	void measure(opercode op, const std::vector<state_idx> &) {
		if(!op.diagonal()) {
			const auto &bond = lat_.bonds[op.bond()];
			const auto &bi = lat_.get_uc_site(bond.i).basis;
			const auto &bj = lat_.get_uc_site(bond.j).basis;

			double j20 = bi.states[op.leg_state(2)].j-bi.states[op.leg_state(0)].j;
			double j31 = bj.states[op.leg_state(3)].j-bj.states[op.leg_state(1)].j;
			
			tmpj_ += lat_.sites[bond.i].sublattice*j20 + lat_.sites[bond.j].sublattice*j31;
		}

		j_ += tmpj_;
		j2_ += tmpj_ * tmpj_;
		n_++;
	}

	void result(loadl::measurements &measure) {
		std::string p = "Sign";

		double norm = 1. / lat_.sites.size();
		j_ *= norm;
		j2_ *= norm * norm;

		measure.add(p + "J", sign_*j_ / n_);
		measure.add(p + "J2", sign_*j2_ / n_);
	}

	void register_evalables(loadl::evaluator &eval) {
		std::string p = "Sign";
		auto unsign = [](const std::vector<std::vector<double>> &obs) {
			return std::vector<double>{obs[0][0]/obs[1][0]};
		};
		
		eval.evaluate( "J", {p+"J", "Sign"}, unsign);
		eval.evaluate( "JVar", {p+"J2", p+"J", "Sign"}, [](const std::vector<std::vector<double>> &obs) {
			return std::vector<double>{(obs[0][0]-obs[1][0]*obs[1][0]/obs[2][0])/obs[2][0]};
		});
	}
};
