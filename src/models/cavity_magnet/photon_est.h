#pragma once

#include "models/cavity_magnet/cavity_magnet.h"
#include "opercode.h"
#include "sse_data.h"
#include <loadleveller/evalable.h>
#include <loadleveller/measurements.h>

#include <iostream>
class photon_est {
private:
	double n_{1};
	std::vector<double> tmpnum_;

	std::vector<double> num_;
	std::vector<double> num2_;

	const double sign_{};
	const cavity_magnet &model_;

public:
	photon_est(const cavity_magnet &model, double sign)
	    : tmpnum_(model.modes.size()), num_(model.modes.size()),
	      num2_(model.modes.size()), sign_{sign}, model_{model} {}

	void init(const std::vector<state_idx> &state) {
		std::fill(tmpnum_.begin(), tmpnum_.end(), 0);

		for(int m = 0; m < static_cast<int>(model_.modes.size()); m++) {
			tmpnum_[m] = model_.get_basis(m).n(state[m]);
			num_[m] = tmpnum_[m];
			num2_[m] = tmpnum_[m] * tmpnum_[m];
		}
	}

	void measure(opercode op, const std::vector<state_idx> &, const sse_data &data) {
		if(!op.diagonal()) {
			const auto &bond = data.get_bond(op.bond());
			int leg_count = data.get_vertex_data(op.bond()).leg_count;

			const auto &leg_state = data.get_vertex_data(op.bond()).get_legstate(op.vertex());

			// photon spin coupling. assumed to be the only offdiagonal photon operator now.
			if(leg_count == 6) {
				const auto &b = model_.get_basis(bond[0]);
				tmpnum_[bond[0]] += b.n(leg_state[leg_count / 2]) - b.n(leg_state[0]);
			}
		}

		for(int m = 0; m < static_cast<int>(model_.modes.size()); m++) {
			num_[m] += tmpnum_[m];
			num2_[m] += tmpnum_[m] * tmpnum_[m];
			n_++;
		}
	}

	void result(loadl::measurements &measure) {
		std::string p = "Sign";

		for(int m = 0; m < static_cast<int>(model_.modes.size()); m++) {
			num_[m] *= sign_ / n_;
			num2_[m] *= sign_ / n_;
		}

		measure.add(p + "PhotonNum", num_);
		measure.add(p + "PhotonNum2", num2_);
	}

	void register_evalables(loadl::evaluator &eval) {
		std::string p = "Sign";

		auto unsign = [](const std::vector<std::vector<double>> &obs) {
			std::vector<double> res(obs[0].size());
			int i{};
			for(auto &r : res) {
				r = obs[0][i] / obs[1][i];
				i++;
			}
			return res;
		};

		for(const auto &obs : {"PhotonNum", "PhotonNum2"}) {
			eval.evaluate(obs, {p + obs, "Sign"}, unsign);
		}

		eval.evaluate("PhotonNumVar", {"SignPhotonNum2", "SignPhotonNum", "Sign"},
		              [](const std::vector<std::vector<double>> &obs) {
			              std::vector<double> res(obs[0].size());
			              int i{};
			              for(auto &r : res) {
				              r = (obs[0][i] - obs[1][i] * obs[1][i]) / obs[2][0];
				              i++;
			              }
			              return res;
		              });
	}
};