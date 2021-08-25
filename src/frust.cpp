#include "frust.h"
#include <algorithm>
#include "vertex_data.h"
#include "latticedef.h"
#include "mag_est.h"
#include "j_est.h"

#include "template_selector.h"

frust::frust(const loadl::parser &p)
	: loadl::mc(p),
	  lat_{lattice_from_param(p, true)}, settings_{p} {
	T_ = param.get<double>("T");
	v_first_.resize(lat_.sites.size());
	v_last_.resize(lat_.sites.size());

}

void frust::print_operators() {
	int p = -1;
	for(auto op : operators_) {
		p++;
		if(op.identity()) {
			continue;
		}
		const auto &b = lat_.bonds[op.bond()];
		const auto &ls = lat_.get_vertex_data(op.bond()).get_legstate(op.vertex());
		const auto &bi = lat_.get_uc_site(b.i).basis.states;
		const auto &bj = lat_.get_uc_site(b.j).basis.states;
		std::cout << fmt::format("{} {}-{}: {}{}->{}{}\n", 4*p, b.i, b.j, bi[ls[0]].name, bj[ls[1]].name, bi[ls[2]].name, bj[ls[3]].name);
	}
	int idx{};
	for(auto s : spin_) {
		std::cout << lat_.get_uc_site(idx).basis.states[s].name << ", ";
	}
	std::cout << "\n";
}

void frust::init() {
	spin_.resize(lat_.sites.size());
	for(size_t i = 0; i < spin_.size(); i++) {
		spin_[i] = lat_.get_uc_site(i).basis.size()*random01();
	}

	operators_.resize(param.get("init_opstring_cutoff",static_cast<int>(lat_.spinhalf_count*T_)));

	const int warmup = 5;
	for(int i = 0; i < warmup; i++) {
		diagonal_update();
	}
}

int frust::worm_traverse() {
	if(noper_ == 0) {
		return 0;
	}
	
	int wormlength{};

	uint32_t v0{};
	do {
		v0 = vertices_.size()*random01();
	} while(vertices_[v0] < 0);

	auto op0 = operators_[v0/4];
	const auto &bond0 = lat_.bonds[op0.bond()];
	int site0 = v0&1 ? bond0.j : bond0.i;
	int wormfunc0 = random01()*lat_.get_uc_site(site0).basis.worms.size();

	uint32_t v = v0;
	int wormfunc = wormfunc0;

	do {
		wormlength++;

		auto &op = operators_[v/4];
		assert(!op.vertex().invalid());
		int leg_in = v%4;
		const auto [leg_out, wormfunc_out, new_vertex] = lat_.get_vertex_data(op.bond()).scatter(op.vertex(), leg_in, wormfunc, random01());
		
		op = opercode{op.bond(), new_vertex};
		
		uint32_t vstep = 4*(v/4)+leg_out;
		const auto &bond = lat_.bonds[op.bond()];
		const auto &site_out = lat_.get_uc_site(leg_out&1 ? bond.j : bond.i);
		if(vstep == v0 && wormfunc_out == site_out.basis.worms[wormfunc0].inverse_idx) {
			break;
		}
		wormfunc = wormfunc_out;
		v = vertices_[vstep];
		assert(vertices_[vstep] != -1);
	} while(v != v0 || wormfunc != wormfunc0);

	return wormlength;
}

std::optional<uint32_t> frust::find_worm_measure_start(int site0, uint32_t &p0, int direction0) const {
	int opsize = operators_.size();

	if(v_first_[site0] < 0) {
		return std::nullopt;
	}

	for(int l = 0; l <= opsize; l++) {
		int p = p0 + direction0*l; 
		if(p < 0) {
			p += opsize;
		}
		if(p >= opsize) {
			p -= opsize;
		}

		auto op = operators_[p];
		if(!op.identity() && (l != 0 || direction0 < 0)) {
			const auto &bond = lat_.bonds[op.bond()];
			if(bond.i == site0 || bond.j == site0) {
				/*if(direction0 == 1) {
					p0 = p0 ? p0-1 : opsize-1;
				}*/
				return 4*p + 2*(direction0<0) + (bond.j==site0);
			}
		}
	}
	assert(false);

	return std::nullopt;
}

void frust::worm_traverse_measure(double &sign, std::vector<double> &corr) {
	const auto matelem_idx = [](bool switcheroo, const site_basis::state &sbefore, const site_basis::state &safter) -> double {
		const auto &sup = switcheroo ? sbefore : safter;
		const auto &sdown = switcheroo ? safter : sbefore;

		if(sdown.j == sup.j && sdown.m == sup.m) {
			return sup.jdim;
		}
		return -1;
	};

	assert(lat_.uc.sites.size() == 1); // not implemented
	
	
	if(noper_ == 0) {
		return;
	}

	uint32_t p0 = operators_.size()*random01();
	int site0 = lat_.sites.size()*random01();
	int direction0 = 1-2*(random01()>0.5);

	auto v0opt = find_worm_measure_start(site0, p0, direction0);
	if(!v0opt) {
		return;
	}
	uint32_t v0 = *v0opt;
	uint32_t v1 = vertices_[v0];
	auto [uc0, x0, y0] = lat_.split_idx(site0);

	const auto &basis0 = lat_.get_uc_site(site0).basis;
	int wormfunc0 = random01()*basis0.worms.size();

	uint32_t v = v0;
	int wormfunc = wormfunc0;

	do {
		auto &op = operators_[v/4];
		assert(!op.vertex().invalid());
		opercode old_op = op;

		int leg_in = v%4;
		const auto &vd = lat_.get_vertex_data(op.bond());
		const auto [leg_out, wormfunc_out, new_vertex] = vd.scatter(op.vertex(), leg_in, wormfunc, random01());
		
		op = opercode{op.bond(), new_vertex};
		
		uint32_t vstep = 4*(v/4)+leg_out;
		const auto &bond = lat_.bonds[op.bond()];
		int site_idx = leg_out&1 ? bond.j : bond.i;
		const auto &site_out = lat_.get_uc_site(site_idx);

		sign *= vd.get_sign(old_op.vertex())*vd.get_sign(op.vertex());

		if(vstep == v0 && wormfunc_out == site_out.basis.worms[wormfunc0].inverse_idx) {
			break;
		}
	
		wormfunc = wormfunc_out;
		uint32_t vnext = vertices_[vstep];
		
		assert(vnext != -1);

		bool up = leg_out > 1;
		if((up  &&     v/4 <= p0 && p0 < vnext/4) ||
		   (!up && vnext/4 <= p0 && p0 < v/4) ||
		   (up  && vnext/4 <=     v/4 && (p0 >=     v/4 || p0 < vnext/4)) ||
		   (!up &&     v/4 <= vnext/4 && (p0 >= vnext/4 || p0 <     v/4))) {
			int state_after_idx = vd.get_legstate(op.vertex())[leg_out];
			int state_before_idx = site_out.basis.worms[site_out.basis.worms[wormfunc].inverse_idx].action[state_after_idx];
			const auto &state_before = site_out.basis.states[state_before_idx];
			const auto &state_after = site_out.basis.states[state_after_idx];	
			int matidx = matelem_idx(!up, state_before, state_after);
			
			if(matidx >= 0 && site_idx != site0) {
				auto op0 = operators_[v0/4];
				auto op1 = operators_[v1/4];
				const auto &ls0 = lat_.get_vertex_data(op0.bond()).get_legstate(op0.vertex());
				const auto &ls1 = lat_.get_vertex_data(op1.bond()).get_legstate(op1.vertex());
				const auto &state_before0 = basis0.states[ls0[v0%4]];
				const auto &state_after0 = basis0.states[ls1[v1%4]];
				double matidx0 = matelem_idx(direction0 > 0, state_before0, state_after0);

				if(matidx0 >= 0) {
					auto [uc, x, y] = lat_.split_idx(site_idx);

					int idx = ((y-y0+lat_.Ly)%lat_.Ly)*lat_.Lx + (x-x0+lat_.Lx)%lat_.Lx;
					if(settings_.loopcorr_as_strucfac) {
						idx = 0;
					}
					corr[4*idx + 2*matidx0 + matidx] += basis0.worms.size()*sign;
				}
			}
		}
		
		v = vnext;
	} while(v != v0 || wormfunc != wormfunc0);
}		
		

void frust::worm_update() {
	if(!is_thermalized() || !settings_.measure_chirality) {
		double wormlength = 1;
		for(int i = 0; i < nworm_; i++) {
			wormlength += worm_traverse();
		}
		
		if(is_thermalized() && noper_ != 0) {
			measure.add("wormlength_fraction", wormlength/noper_);
		}
		
		wormlength /= ceil(nworm_);

		if(!is_thermalized()) {
			avgwormlen_ += 0.01*(wormlength-avgwormlen_);
			double target_worms = param.get<double>("wormlength_fraction", 2.0)*noper_/avgwormlen_;
			nworm_ += 0.01*(target_worms-nworm_)+tanh(target_worms-nworm_);
			nworm_ = std::clamp(nworm_, 1., 1.+noper_/2.);
		}
	} else {
		int size = settings_.loopcorr_as_strucfac ? 1 : lat_.sites.size();
		std::vector<double> corr(4 * size);
		double sign = measure_sign();

		for(int i = 0; i < nworm_; i++) {
			worm_traverse_measure(sign, corr);
		}	
		for(auto &c : corr) {
			c *= 1./ceil(nworm_);
		}

		measure.add("SignJDimOffCorr", corr);
	}

	for(size_t i = 0; i < spin_.size(); i++) {
		if(v_first_[i] < 0) {
			spin_[i] = lat_.get_uc_site(i).basis.size() * random01();
		} else {
			uint32_t v = v_first_[i];
			auto op = operators_[v/4];
			int leg = v%4;
			const auto &leg_state = lat_.get_vertex_data(op.bond()).get_legstate(op.vertex());
			spin_[i] = leg_state[leg];
		}
	}

}

void frust::make_vertex_list() {
	vertices_.resize(operators_.size() * 4);
	std::fill(vertices_.begin(), vertices_.end(), -1);
	std::fill(v_first_.begin(), v_first_.end(), -1);
	std::fill(v_last_.begin(), v_last_.end(), -1);

	for(size_t p = 0; p < operators_.size(); p++) {
		auto op = operators_[p];
		if(op.identity()) {
			continue;
		}
		int v0 = 4 * p;
		const auto &b = lat_.bonds[op.bond()];
		int v1 = v_last_[b.i];
		int v2 = v_last_[b.j];
		if(v1 != -1) {
			vertices_[v1] = v0;
			vertices_[v0] = v1;
		} else {
			v_first_[b.i] = v0;
		}
		if(v2 != -1) {
			vertices_[v2] = v0 + 1;
			vertices_[v0 + 1] = v2;
		} else {
			v_first_[b.j] = v0 + 1;
		}
		v_last_[b.i] = v0 + 2;
		v_last_[b.j] = v0 + 3;
	}

	for(size_t i = 0; i < v_first_.size(); i++) {
		int f = v_first_[i];
		if(f != -1) {
			int l = v_last_[i];
			vertices_[f] = l;
			vertices_[l] = f;
		}
	}

}

void frust::diagonal_update() {
	if(noper_ >= operators_.size() * 0.5) {
		if(is_thermalized()) {
			std::cout << "Warning: spin array resized after thermalization\n";
		}
		operators_.resize(operators_.size() * 1.5 + 100, opercode::make_identity());
	}

	auto tmpspin = spin_;

	double opersize = operators_.size();
	double p_make_bond_raw = lat_.bonds.size()*1./T_;
	double p_remove_bond_raw = T_/lat_.bonds.size();

	for(auto &op : operators_) {

		if(op.identity()) {
			uint32_t bond = random01() * lat_.bonds.size();
			const auto &b = lat_.bonds[bond];
			state_idx si = spin_[b.i];
			state_idx sj = spin_[b.j];

			vertexcode newvert = lat_.get_vertex_data(bond).get_diagonal_vertex(si, sj);
			double weight = lat_.get_vertex_data(bond).get_weight(newvert);
			double p_make_bond = p_make_bond_raw/(opersize-noper_);

			if(random01() < p_make_bond * weight) {
				op = opercode{bond, newvert};
				noper_++;
			}
		} else {
			int bond = op.bond();
			const auto &vertdata = lat_.get_vertex_data(bond);
			if(op.diagonal()) {
				double weight = vertdata.get_weight(op.vertex());
				double p_remove_bond = (opersize-noper_+1)*p_remove_bond_raw;
				if(random01()*weight < p_remove_bond) {
					op = opercode::make_identity();
					noper_--;
				}
			} else {
				const auto &b = lat_.bonds[bond];
				const auto &leg_state = vertdata.get_legstate(op.vertex());
				spin_[b.i] = leg_state[2];
				spin_[b.j] = leg_state[3];
			}
		}
	}
	assert(tmpspin == spin_);
}

void frust::do_update() {
	diagonal_update();
	make_vertex_list();
	worm_update();
}

template<typename... Ests>
void frust::opstring_measurement(Ests... ests) {
	auto x1 = {(ests.init(spin_), 0)...};
	(void)x1;

	int64_t n = 0;
	for(auto op : operators_) {
		if(op.identity()) {
			continue;
		}

		if(!op.diagonal()) {
			const auto &bond = lat_.bonds[op.bond()];
			const auto &leg_state = lat_.get_vertex_data(op.bond()).get_legstate(op.vertex());
			spin_[bond.i] = leg_state[2];
			spin_[bond.j] = leg_state[3];
		}

		if(n < noper_) {
			auto x2 = {(ests.measure(op, spin_), 0)...};
			(void)x2;
		}
		n++;
	}

	auto x3 = {(ests.result(measure), 0)...};
	(void)x3;
}

double frust::measure_sign() const {
	int sign{};

	for(auto op : operators_) {
		if(!op.identity()) {
			sign += lat_.get_vertex_data(op.bond()).get_sign(op.vertex()) < 0;
		}
	}
	return (sign&1) ? -1 : 1;
}

void frust::do_measurement() {
	double sign = measure_sign();

	auto obs = std::tuple{
	    j_est{lat_, sign, settings_.measure_jcorrlen},
	    mag_est<1, 1>{lat_, T_, sign},
	    mag_est<1, -1>{lat_, T_, sign},
	    mag_est<-1, 1>{lat_, T_, sign},
	    mag_est<-1, -1>{lat_, T_, sign},
	};

	std::array<bool, 5> flags = {
	    settings_.measure_j || settings_.measure_chirality,     settings_.measure_mag,     settings_.measure_sxmag,
	    settings_.measure_symag, settings_.measure_sxsymag,
	};

	template_select([&](auto... vals) { opstring_measurement(vals...); }, obs, flags);

	measure.add("Sign", sign);
	measure.add("nOper", static_cast<double>(noper_));
	measure.add("SignNOper", sign*static_cast<double>(noper_));
	measure.add("SignNOper2", sign*static_cast<double>(noper_)*static_cast<double>(noper_));

	measure.add("SignEnergy", sign*(-static_cast<double>(noper_) * T_ - lat_.energy_offset) / lat_.spinhalf_count);
}

void frust::checkpoint_write(const loadl::iodump::group &out) {
	std::vector<unsigned int> saveops(operators_.size());
	std::transform(operators_.begin(), operators_.end(), saveops.begin(), [](opercode op) {
		return op.code();
	});
	
	out.write("noper", noper_);
	out.write("avgwormlen", avgwormlen_);
	out.write("nworm", nworm_);
	out.write("operators", saveops);
	out.write("spin", spin_);
	out.write("version", dump_version_);
}

void frust::checkpoint_read(const loadl::iodump::group &in) {
	std::vector<unsigned int> saveops;

	in.read("noper", noper_);
	in.read("avgwormlen", avgwormlen_);
	in.read("nworm", nworm_);
	in.read("operators", saveops);

	operators_.resize(saveops.size());
	std::transform(saveops.begin(), saveops.end(), operators_.begin(), [](uint32_t c) { return opercode{c}; });
	
	in.read("spin", spin_);

	int version;
	in.read("version", version);
	if(version != dump_version_) {
		throw std::runtime_error("dump file does not fit program version");
	}
}

template<int CoeffDim0, int CoeffDim1, int CoeffM>
static double corrfunc_matrix(const std::vector<double> &corr, int idx) {
	double result = 0;
	for(int jdimi = 0; jdimi < 2; jdimi++) {
		for(int jdimj = 0; jdimj < 2; jdimj++) {
			int matidx = 2*jdimi + jdimj;
			double coeff = (CoeffDim0 + (CoeffDim1-CoeffDim0) * jdimi) * (CoeffDim0 + (CoeffDim1-CoeffDim0) * jdimj);
			result += coeff * corr[4*idx + matidx];
		}
	}
	return result;
}

void frust::register_evalables(loadl::evaluator &eval, const loadl::parser &p) {
	auto unsign = [](const std::vector<std::vector<double>> &obs) {
		std::vector<double> result = obs[0];
		for(auto &r : result) {
			r /= obs[1][0];
		}
		
		return result;
	};

	double T = p.get<double>("T");
	measurement_settings settings{p};
	lattice lat = lattice_from_param(p, false);

	if(settings.measure_j || settings.measure_chirality) {
		j_est{lat,0,settings.measure_jcorrlen}.register_evalables(eval);
	}
	
	if(settings.measure_mag) {
		mag_est<1,1>{lat, T, 0}.register_evalables(eval);
	}

	if(settings.measure_sxmag) {
		mag_est<-1,1>{lat, T, 0}.register_evalables(eval);
	}
	
	if(settings.measure_symag) {
		mag_est<1,-1>{lat, T, 0}.register_evalables(eval);
	}

	if(settings.measure_sxsymag) {
		mag_est<-1,-1>{lat, T, 0}.register_evalables(eval);
	}

	if(settings.measure_chirality) {
		if(settings.loopcorr_as_strucfac) {
			eval.evaluate("NematicityStruc", {"SignJDimOffCorr", "SignChiralityOnsite", "Sign", "SignNematicityDiagStruc"}, [](const std::vector<std::vector<double>> &obs) {
				double result = 9.0/16.0*(obs[1][0]+corrfunc_matrix<1,1,1>(obs[0], 0) + obs[3][0])/obs[2][0];
				return std::vector<double>{result};
			});
			eval.evaluate("TauZStruc", {"SignJDimOffCorr", "SignChiralityOnsite", "Sign"}, [](const std::vector<std::vector<double>> &obs) {
				double result = (obs[1][0] + corrfunc_matrix<-1,1,1>(obs[0], 0))/obs[2][0];
				return std::vector<double>{result};
			});
			eval.evaluate("TauYStruc", {"SignJDimOffCorr", "SignChiralityOnsite", "Sign"}, [](const std::vector<std::vector<double>> &obs) {
				double result = (obs[1][0] + corrfunc_matrix<1,1,1>(obs[0], 0))/obs[2][0];

				return std::vector<double>{result};
			});
		} else {
			eval.evaluate("NematicityOffCorr", {"SignJDimOffCorr", "SignChiralityOnsite", "Sign"}, [](const std::vector<std::vector<double>> &obs) {
				std::vector<double> result(obs[0].size()/4, 0);
				for(int i = 1; i < result.size(); i++) {
					result[i] = 9.0/16.0 * corrfunc_matrix<1,1,1>(obs[0],i)/obs[2][0];
				}
				result[0] += 9.0/16.0 * obs[1][0]/obs[2][0];
				return result;
			});
			eval.evaluate("TauZ", {"SignJDimOffCorr", "SignChiralityOnsite", "Sign"}, [](const std::vector<std::vector<double>> &obs) {
				std::vector<double> result(obs[0].size()/4, 0);
				for(int i = 1; i < result.size(); i++) {
					result[i] = corrfunc_matrix<-1,1,1>(obs[0],i)/obs[2][0];
				}
				result[0] += obs[1][0]/obs[2][0];
				return result;
			});
			eval.evaluate("TauY", {"SignJDimOffCorr", "SignChiralityOnsite", "Sign"}, [](const std::vector<std::vector<double>> &obs) {
				std::vector<double> result(obs[0].size()/4, 0);
				for(int i = 1; i < result.size(); i++) {
					result[i] = corrfunc_matrix<1,1,1>(obs[0], i)/obs[2][0];
				}
				result[0] += obs[1][0]/obs[2][0];

				return result;
			});

		}
	}
	
	eval.evaluate("Energy", {"SignEnergy", "Sign"}, unsign);
	eval.evaluate("SpecificHeat", {"SignNOper2", "SignNOper", "Sign"}, [&](const std::vector<std::vector<double>> &obs) {
		double sn2 = obs[0][0];
		double sn = obs[1][0];
		double sign = obs[2][0];

		return std::vector<double>{(sn2/sign - sn*sn/sign/sign - sn/sign)/lat.spinhalf_count};
	});
}

double frust::pt_weight_ratio(const std::string &param_name, double new_param) {
	// we only implement temperature swaps
	if(param_name != "T") {
		return 0;
	}

	return -noper_*(log(new_param/T_));
}

void frust::pt_update_param(const std::string &param_name, double new_param) {
	if(param_name == "T") {
		T_ = new_param;
	}
}
