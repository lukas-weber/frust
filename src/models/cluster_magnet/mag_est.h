#pragma once

#include "cluster_magnet.h"
#include "opercode.h"
#include "sse_data.h"
#include <loadleveller/evalable.h>
#include <loadleveller/measurements.h>

template<int SignX, int SignY, int SignUC>
class mag_est {
private:
	double n_{1};

	double tmpmag_{};
	double mag_{};
	double absmag_{};
	double mag2_{};
	double mag4_{};

	const cluster_magnet &model_;
	double T_{};
	double sign_{};

	double stag_sign(uint32_t idx) const {
		double sign = 1;

		if(SignUC < 0) {
			sign *= model_.lat.site_sublattice_sign(idx);
		}

		idx /= model_.lat.uc.sites.size();
		if(SignX < 0) {
			sign *= 2. * ((idx % model_.lat.Lx) & 1) - 1.;
		}
		if(SignY < 0) {
			sign *= 2. * ((idx / model_.lat.Lx) & 1) - 1.;
		}

		return sign;
	}

	const std::string prefix() const {
		std::string p;
		if(SignX < 0) {
			p += "StagX";
		}
		if(SignY < 0) {
			p += "StagY";
		}
		if(SignUC < 0) {
			p += "StagUC";
		}
		return p;
	}

public:
	mag_est(const model &model, double T, double sign)
	    : model_{static_cast<const cluster_magnet &>(model)}, T_{T}, sign_{sign} {}

	void init(const std::vector<state_idx> &spin) {
		tmpmag_ = 0;

		for(uint32_t i = 0; i < spin.size(); i++) {
			tmpmag_ += stag_sign(i) * model_.get_site(i).basis.m(spin[i]);
		}

		mag_ = tmpmag_;
		absmag_ = fabs(mag_);
		mag2_ = mag_ * mag_;
		mag4_ = mag2_ * mag2_;
	}

	void measure(opercode op, const std::vector<state_idx> &, const sse_data &data) {
		if(SignX == 1 && SignY == 1) {
			return;
		}

		if(!op.diagonal()) {
			const auto &bond = model_.lat.bonds[op.bond()];
			const auto &bi = model_.get_site(bond.i).basis;
			const auto &bj = model_.get_site(bond.j).basis;

			const auto &leg_state = data.get_vertex_data(op.bond()).get_legstate(op.vertex());

			double m20 = bi.m(leg_state[2]) - bi.m(leg_state[0]);
			double m31 = bj.m(leg_state[3]) - bj.m(leg_state[1]);

			tmpmag_ += stag_sign(bond.i) * m20 + stag_sign(bond.j) * m31;
		}

		mag_ += tmpmag_;
		absmag_ += fabs(tmpmag_);
		mag2_ += tmpmag_ * tmpmag_;
		mag4_ += tmpmag_ * tmpmag_ * tmpmag_ * tmpmag_;
		n_++;
	}

	void result(loadl::measurements &measure) {
		std::string p = "Sign";
		p += prefix();

		double norm = 1. / model_.spinhalf_count;
		mag_ *= norm;
		absmag_ *= norm;
		mag2_ *= norm * norm;
		mag4_ *= norm * norm * norm * norm;

		measure.add(p + "Mag", sign_ * mag_ / n_);
		measure.add(p + "AbsMag", sign_ * absmag_ / n_);
		measure.add(p + "Mag2", sign_ * mag2_ / n_);
		measure.add(p + "Mag4", sign_ * mag4_ / n_);

		double chi = 1 / T_ / (n_ + 1) / n_ * (mag_ * mag_ + mag2_) * model_.spinhalf_count;
		measure.add(p + "MagChi", sign_ * chi);
	}

	void register_evalables(loadl::evaluator &eval) {
		std::string p = "Sign";
		std::string p2 = prefix();

		auto unsign = [](const std::vector<std::vector<double>> &obs) {
			return std::vector<double>{obs[0][0] / obs[1][0]};
		};

		eval.evaluate(p2 + "Mag", {p + p2 + "Mag", "Sign"}, unsign);
		eval.evaluate(p2 + "Mag2", {p + p2 + "Mag2", "Sign"}, unsign);
		eval.evaluate(p2 + "Mag4", {p + p2 + "Mag4", "Sign"}, unsign);

		eval.evaluate(p2 + "BinderRatio", {p + p2 + "Mag4", p + p2 + "Mag2", "Sign"},
		              [](const std::vector<std::vector<double>> &obs) {
			              double smag4 = obs[0][0];
			              double smag2 = obs[1][0];
			              double sign = obs[2][0];

			              return std::vector<double>{smag2 * smag2 / smag4 / sign};
		              });

		eval.evaluate(p2 + "MagChi", {p + p2 + "MagChi", "Sign"}, unsign);
	}
};
