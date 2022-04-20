#pragma once

#include "cluster_magnet.h"
#include "opercode.h"
#include "sse_data.h"
#include <loadleveller/evalable.h>
#include <loadleveller/measurements.h>

namespace mag_sign {
constexpr uint32_t none = 0;
constexpr uint32_t x = 1;
constexpr uint32_t y = 2;
constexpr uint32_t uc = 4;
};

template<uint32_t Signs, typename Model>
class mag_est {
private:
	double n_{1};

	double tmpmag_{};
	double mag_{};
	double absmag_{};
	double mag2_{};
	double mag4_{};

	const Model &model_;
	double T_{};
	double sign_{};

	double stag_sign(lattice::site_idx idx) const {
		double sign = 1;

		if(Signs & mag_sign::uc) {
			sign *= model_.lat.site_sublattice_sign(idx);
		}

		idx /= model_.lat.uc.sites.size();
		if(Signs & mag_sign::x) {
			sign *= 2. * ((idx % model_.lat.Lx) & 1) - 1.;
		}
		if(Signs & mag_sign::y) {
			sign *= 2. * ((idx / model_.lat.Lx) & 1) - 1.;
		}

		return sign;
	}

	const std::string prefix() const {
		std::string p;
		if(Signs & mag_sign::x) {
			p += "StagX";
		}
		if(Signs & mag_sign::y) {
			p += "StagY";
		}
		if(Signs & mag_sign::uc) {
			p += "StagUC";
		}
		return p;
	}

public:
	mag_est(const Model &model, double T, double sign) : model_{model}, T_{T}, sign_{sign} {}

	void init(const std::vector<state_idx> &spin) {
		tmpmag_ = 0;

		for(int i = 0; i < static_cast<int>(spin.size()); i++) {
			if(auto site = model_.get_lattice_site_idx(i)) {
				tmpmag_ += stag_sign(*site) * model_.get_basis(i).m(spin[i]);
			}
		}

		mag_ = tmpmag_;
		absmag_ = fabs(mag_);
		mag2_ = mag_ * mag_;
		mag4_ = mag2_ * mag2_;
	}

	void measure(opercode op, const std::vector<state_idx> &, const sse_data &data) {
		if(!(Signs & mag_sign::x) && !(Signs & mag_sign::y)) {
			return;
		}

		if(!op.diagonal()) {
			const auto &bond = data.get_bond(op.bond());

			const auto &leg_state = data.get_vertex_data(op.bond()).get_legstate(op.vertex());

			int leg_count = data.get_vertex_data(op.bond()).leg_count;
			double mdiff = 0;
			for(int l = 0; l < leg_count / 2; l++) {
				const auto &b = model_.get_basis(bond[l]);
				if(auto site = model_.get_lattice_site_idx(bond[l])) {
					mdiff += stag_sign(*site) * b.m(leg_state[leg_count / 2 + l]) - b.m(l);
				}
			}
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

		double norm = 1. / model_.normalization_site_count();
		mag_ *= norm;
		absmag_ *= norm;
		mag2_ *= norm * norm;
		mag4_ *= norm * norm * norm * norm;

		measure.add(p + "Mag", sign_ * mag_ / n_);
		measure.add(p + "AbsMag", sign_ * absmag_ / n_);
		measure.add(p + "Mag2", sign_ * mag2_ / n_);
		measure.add(p + "Mag4", sign_ * mag4_ / n_);

		double chi =
		    1 / T_ / (n_ + 1) / n_ * (mag_ * mag_ + mag2_) * model_.normalization_site_count();
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
