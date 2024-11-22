#include "cluster_magnet.h"
#include "../common/mag_est.h"
#include "j_est.h"

#include "util/kronecker_product.h"
#include <Eigen/Dense>

static Eigen::MatrixXd onsite_term(const cluster_site &s) {
	const auto &b = s.basis;

	Eigen::MatrixXd res = res.Zero(b.size(), b.size());

	if(static_cast<int>(s.Jin.size()) != b.nspinhalfs * (b.nspinhalfs - 1) / 2) {
		throw std::runtime_error{"Jin does not have the right shape"};
	}

	int idx = 0;
	for(int i = 0; i < b.nspinhalfs; i++) {
		auto spini = b.spinop(i);
		res += spini[1] * s.h;
		for(int j = 0; j < i; j++) {
			auto spinj = b.spinop(j);
			res += s.Jin[idx] *
			       (0.5 * (spini[0] * spinj[0].transpose() + spini[0].transpose() * spinj[0]) +
			        spini[1] * spinj[1]);
			idx++;
		}
	}
	return res;
}

static Eigen::MatrixXd bond_term(const cluster_bond &b, const cluster_site &si,
                                 const cluster_site &sj) {
	assert(static_cast<int>(b.J.size()) == si.basis.nspinhalfs * sj.basis.nspinhalfs);

	int dim = si.basis.size() * sj.basis.size();
	Eigen::MatrixXd res = Eigen::MatrixXd::Zero(dim, dim);

	for(int i = 0; i < si.basis.nspinhalfs; i++) {
		auto spini = si.basis.spinop(i);
		for(int j = 0; j < sj.basis.nspinhalfs; j++) {
			auto spinj = sj.basis.spinop(j);
			res += b.J[sj.basis.nspinhalfs * i + j] * scalar_product(spini, spinj);
		}
	}

	return res;
}

cluster_magnet::cluster_magnet(const lattice &lat, const std::vector<cluster_site> &sites,
                               const std::vector<cluster_bond> &bonds,
                               const cluster_magnet_measurement_settings &settings)
    : model{model::model_type::cluster_magnet}, lat{lat}, settings{settings}, bonds_{bonds},
      sites_{sites} {
	assert(bonds_.size() == lat.uc.bonds.size());
	assert(sites_.size() == lat.uc.sites.size());

	for(const auto &site : sites_) {
		spinhalf_count += site.basis.nspinhalfs;
	}

	spinhalf_count *= lat.Lx * lat.Ly;
}

sse_data cluster_magnet::generate_sse_data() const {
	std::vector<vertex_data> vert_data;
	for(int b = 0; b < static_cast<int>(bonds_.size()); b++) {
		const auto &bond = lat.uc.bonds[b];
		const auto &si = sites_[bond.i];
		const auto &sj = sites_[bond.j.uc];
		int dim_i = si.basis.size();
		int dim_j = sj.basis.size();

		int dim = dim_i * dim_j;
		Eigen::MatrixXd H = Eigen::MatrixXd::Zero(dim, dim);
		auto Idi = Eigen::MatrixXd::Identity(dim_i, dim_i);
		auto Idj = Eigen::MatrixXd::Identity(dim_j, dim_j);
		H += kronecker_prod(onsite_term(si), Idj) / lat.uc.sites[bond.i].coordination;
		H += kronecker_prod(Idi, onsite_term(sj)) / lat.uc.sites[bond.j.uc].coordination;
		H += bond_term(bonds_[b], si, sj);

		vert_data.push_back(vertex_data{{dim_i, dim_j}, H});
	}
	std::vector<sse_data::site> sites;
	std::vector<sse_data::bond> bonds;
	for(int i = 0; i < lat.Lx * lat.Ly; i++) {
		std::transform(sites_.begin(), sites_.end(), std::back_inserter(sites),
		               [](const auto &s) { return sse_data::site{s.basis.size()}; });
	}
	std::transform(lat.bonds.begin(), lat.bonds.end(), std::back_inserter(bonds),
	               [](const auto &b) {
		               return sse_data::bond{b.type, {b.i, b.j}};
	               });

	return sse_data{vert_data, sites, bonds};
}

void cluster_magnet::to_json(nlohmann::json &out) const {
	out["model"] = "cluster_magnet";

	lat.to_json(out);

	for(int i = 0; i < lat.Lx * lat.Ly * static_cast<int>(lat.uc.sites.size()); i++) {
		auto &site = out["sites"][i];
		const auto &s = get_site(i);
		site["nspinhalfs"] = s.basis.nspinhalfs;
		site["Jin"] = s.Jin;
		site["h"] = s.h;
	}

	for(int i = 0; i < static_cast<int>(lat.bonds.size()); i++) {
		out["bonds"][i]["J"] = get_bond(i).J;
	}
}

template<int CoeffDim0, int CoeffDim1, int CoeffM>
static double corrfunc_matrix(const std::vector<double> &corr, int idx) {
	double result = 0;
	for(int jdimi = 0; jdimi < 2; jdimi++) {
		for(int jdimj = 0; jdimj < 2; jdimj++) {
			int matidx = 2 * jdimi + jdimj;
			double coeff = (CoeffDim0 + (CoeffDim1 - CoeffDim0) * jdimi) *
			               (CoeffDim0 + (CoeffDim1 - CoeffDim0) * jdimj);
			result += coeff * corr[4 * idx + matidx];
		}
	}
	return result;
}

void cluster_magnet::register_evalables(loadl::evaluator &eval, double T) const {
	using M = cluster_magnet;

	if(settings.measure_j || settings.measure_chirality) {
		j_est{*this, 0, settings.measure_jcorrlen}.register_evalables(eval);
	}

	if(settings.measure_mag) {
		mag_est<mag_sign::none, M>{*this, T, 0}.register_evalables(eval);
	}

	if(settings.measure_sxmag) {
		mag_est<mag_sign::x, M>{*this, T, 0}.register_evalables(eval);
	}

	if(settings.measure_symag) {
		mag_est<mag_sign::y, M>{*this, T, 0}.register_evalables(eval);
	}

	if(settings.measure_sxsymag) {
		mag_est<mag_sign::x | mag_sign::y, M>{*this, T, 0}.register_evalables(eval);
	}

	if(settings.measure_sxsucmag) {
		mag_est<mag_sign::x | mag_sign::uc, M>{*this, T, 0}.register_evalables(eval);
	}

	if(settings.measure_chirality) {
		if(settings.loopcorr_as_strucfac) {
			eval.evaluate(
			    "NematicityStruc",
			    {"SignJDimOffCorr", "SignChiralityOnsite", "Sign", "SignNematicityDiagStruc"},
			    [](const std::vector<std::vector<double>> &obs) {
				    double result = 9.0 / 16.0 *
				                    (obs[1][0] + corrfunc_matrix<1, 1, 1>(obs[0], 0) + obs[3][0]) /
				                    obs[2][0];
				    return std::vector<double>{result};
			    });
			eval.evaluate("TauZStruc", {"SignJDimOffCorr", "SignChiralityOnsite", "Sign"},
			              [](const std::vector<std::vector<double>> &obs) {
				              double result =
				                  (obs[1][0] + corrfunc_matrix<-1, 1, 1>(obs[0], 0)) / obs[2][0];
				              return std::vector<double>{result};
			              });
			eval.evaluate("TauYStruc", {"SignJDimOffCorr", "SignChiralityOnsite", "Sign"},
			              [](const std::vector<std::vector<double>> &obs) {
				              double result =
				                  (obs[1][0] + corrfunc_matrix<1, 1, 1>(obs[0], 0)) / obs[2][0];

				              return std::vector<double>{result};
			              });
		} else {
			eval.evaluate("NematicityOffCorr", {"SignJDimOffCorr", "SignChiralityOnsite", "Sign"},
			              [](const std::vector<std::vector<double>> &obs) {
				              std::vector<double> result(obs[0].size() / 4, 0);
				              for(int i = 1; i < static_cast<int>(result.size()); i++) {
					              result[i] =
					                  9.0 / 16.0 * corrfunc_matrix<1, 1, 1>(obs[0], i) / obs[2][0];
				              }
				              result[0] += 9.0 / 16.0 * obs[1][0] / obs[2][0];
				              return result;
			              });
			eval.evaluate("TauZ", {"SignJDimOffCorr", "SignChiralityOnsite", "Sign"},
			              [](const std::vector<std::vector<double>> &obs) {
				              std::vector<double> result(obs[0].size() / 4, 0);
				              for(int i = 1; i < static_cast<int>(result.size()); i++) {
					              result[i] = corrfunc_matrix<-1, 1, 1>(obs[0], i) / obs[2][0];
				              }
				              result[0] += obs[1][0] / obs[2][0];
				              return result;
			              });
			eval.evaluate("TauY", {"SignJDimOffCorr", "SignChiralityOnsite", "Sign"},
			              [](const std::vector<std::vector<double>> &obs) {
				              std::vector<double> result(obs[0].size() / 4, 0);
				              for(int i = 1; i < static_cast<int>(result.size()); i++) {
					              result[i] = corrfunc_matrix<1, 1, 1>(obs[0], i) / obs[2][0];
				              }
				              result[0] += obs[1][0] / obs[2][0];

				              return result;
			              });
		}
	}
}
