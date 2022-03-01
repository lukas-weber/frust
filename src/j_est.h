#pragma once

#include "opercode.h"
#include "sse_data.h"
#include "models/model.h"
#include "models/cluster_magnet.h"
#include <loadleveller/evalable.h>
#include <loadleveller/measurements.h>

class j_est {
private:
	double n_{1};

	double tmpj_{};
	double tmpjdim_{};
	double tmpnemdiag_{};
	std::complex<double> tmpjstruc_{};
	std::complex<double> tmpjstruc2_{};
	double j_{};
	double j2_{};
	double jdim_{};
	double nemdiag_{};

	double jstruc_{};
	double jstruc2_{};

	const bool measure_corrlen_{};
	const cluster_magnet &model_;
	double sign_{};
	const double corrq_{};

public:
	j_est(const model &model, double sign, bool measure_corrlen)
	    : measure_corrlen_(measure_corrlen), model_{static_cast<const cluster_magnet &>(model)}, sign_{sign}, corrq_{2 * M_PI / model_.lat.Lx} {}

	void init(const std::vector<state_idx> &spin) {
		using namespace std::complex_literals;
		tmpj_ = 0;
		tmpjdim_ = 0;
		tmpnemdiag_ = 0;
		int uc_size = model_.lat.uc.sites.size();

		int i = 0;
		for(auto s : spin) {
			double j = model_.get_site(i).basis.states[s].j;
			double jdim = model_.get_site(i).basis.states[s].jdim;
			tmpj_ += j;
			tmpjdim_ += jdim;

			tmpnemdiag_ += (2 * jdim - 1) * (1.5 - j);

			if(measure_corrlen_) {
				double xi = (i / uc_size) % model_.lat.Lx;
				tmpjstruc_ += std::exp(1i * corrq_ * xi) * j;
				tmpjstruc2_ += std::exp(2i * corrq_ * xi) * j;
			}
			i++;
		}

		j_ = tmpj_;
		jdim_ = tmpjdim_;
		nemdiag_ = tmpnemdiag_ * tmpnemdiag_;
		j2_ = j_ * j_;

		if(measure_corrlen_) {
			jstruc_ = std::real(tmpjstruc_ * std::conj(tmpjstruc_));
			jstruc2_ = std::real(tmpjstruc2_ * std::conj(tmpjstruc2_));
		}
	}

	void measure(opercode op, const std::vector<state_idx> &, const sse_data &data) {
		using namespace std::complex_literals;
		int uc_size = model_.lat.uc.sites.size();
		if(!op.diagonal()) {
			const auto &bond = model_.lat.bonds[op.bond()];
			const auto &bi = model_.get_site(bond.i).basis;
			const auto &bj = model_.get_site(bond.j).basis;

			const auto &leg_state = data.get_vertex_data(op.bond()).get_legstate(op.vertex());

			double j20 = bi.states[leg_state[2]].j - bi.states[leg_state[0]].j;
			double j31 = bj.states[leg_state[3]].j - bj.states[leg_state[1]].j;

			tmpj_ += j20 + j31;

			double jdim20 = bi.states[leg_state[2]].jdim - bi.states[leg_state[0]].jdim;
			double jdim31 = bj.states[leg_state[3]].jdim - bj.states[leg_state[1]].jdim;

			tmpjdim_ += jdim20 + jdim31;

			double nem20 =
			    (2 * bi.states[leg_state[2]].jdim - 1) * (1.5 - bi.states[leg_state[2]].j) -
			    (2 * bi.states[leg_state[0]].jdim - 1) * (1.5 - bi.states[leg_state[0]].j);
			double nem31 =
			    (2 * bj.states[leg_state[3]].jdim - 1) * (1.5 - bj.states[leg_state[3]].j) -
			    (2 * bj.states[leg_state[1]].jdim - 1) * (1.5 - bj.states[leg_state[1]].j);

			tmpnemdiag_ += nem20 + nem31;

			if(measure_corrlen_) {
				// FIXME: should it not be Ly?
				double xi = (bond.i / uc_size) % model_.lat.Lx;
				double xj = (bond.j / uc_size) % model_.lat.Lx;
				std::complex<double> jq20 = std::exp(1i * corrq_ * xi) *
				                            (bi.states[leg_state[2]].j - bi.states[leg_state[0]].j);
				std::complex<double> jq31 = std::exp(1i * corrq_ * xj) *
				                            (bj.states[leg_state[3]].j - bj.states[leg_state[1]].j);
				tmpjstruc_ += jq20 + jq31;
				tmpjstruc2_ +=
				    std::exp(1i * corrq_ * xi) * jq20 + std::exp(1i * corrq_ * xj) * jq31;
			}
		}

		j_ += tmpj_;
		jdim_ += tmpjdim_;
		nemdiag_ += tmpnemdiag_ * tmpnemdiag_;
		j2_ += tmpj_ * tmpj_;

		if(measure_corrlen_) {
			jstruc_ += std::real(tmpjstruc_ * std::conj(tmpjstruc_));
			jstruc2_ += std::real(tmpjstruc2_ * std::conj(tmpjstruc2_));
		}
		n_++;
	}

	void result(loadl::measurements &measure) {
		std::string p = "Sign";

		double norm = 1. / model_.lat.site_count();
		j_ *= norm;
		jdim_ *= norm;
		nemdiag_ *= norm;
		j2_ *= norm * norm;
		jstruc_ *= norm * norm; // unusual normalization for the structure factor but well…
		jstruc2_ *= norm * norm;

		measure.add(p + "JDim", sign_ * jdim_ / n_);
		measure.add(p + "NematicityDiagStruc", sign_ * nemdiag_ / n_);
		measure.add(p + "J", sign_ * j_ / n_);
		measure.add(p + "J2", sign_ * j2_ / n_);
		if(measure_corrlen_) {
			measure.add(p + "JStruc1", sign_ * jstruc_ / n_);
			measure.add(p + "JStruc2", sign_ * jstruc2_ / n_);
		}

		measure.add(p + "ChiralityOnsite", sign_ * (1.5 - j_ / n_));
	}

	void register_evalables(loadl::evaluator &eval) {
		std::string p = "Sign";
		auto unsign = [](const std::vector<std::vector<double>> &obs) {
			return std::vector<double>{obs[0][0] / obs[1][0]};
		};

		eval.evaluate("JDim", {p + "JDim", "Sign"}, unsign);
		eval.evaluate("NematicityDiagStruc", {p + "NematicityDiagStruc", "Sign"}, unsign);
		eval.evaluate("J", {p + "J", "Sign"}, unsign);
		eval.evaluate("J2", {p + "J2", "Sign"}, unsign);
		if(measure_corrlen_) {
			eval.evaluate("JStruc1", {p + "JStruc1", "Sign"}, unsign);
			eval.evaluate("JStruc2", {p + "JStruc2", "Sign"}, unsign);
			eval.evaluate("JCorrLen", {p + "JStruc1", p + "JStruc2"},
			              [&](const std::vector<std::vector<double>> &obs) {
				              double Sq1 = obs[0][0];
				              double Sq2 = obs[1][0];
				              double r = Sq1 / Sq2;
				              return std::vector<double>{
				                  1 / corrq_ * std::sqrt(std::max((r - 1) / (4 - r), 0.0))};
			              });
		}
		eval.evaluate("JVar", {p + "J2", p + "J", "Sign"},
		              [](const std::vector<std::vector<double>> &obs) {
			              return std::vector<double>{
			                  (obs[0][0] - obs[1][0] * obs[1][0] / obs[2][0]) / obs[2][0]};
		              });
	}
};
