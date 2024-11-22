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

	const int warmup = 5;
	for(int i = 0; i < warmup; i++) {
		diagonal_update();
	}

	measure.register_observable("TauZProb1", 100000);
	measure.register_observable("TauZProb2", 100000);
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

	for(int l = 0; l < opsize; l++) {
		int p = (p0 + direction0*l)%opsize; 

		auto op = operators_[p];
		if(!op.identity()) {
			const auto &bond = lat_.bonds[op.bond()];
			if(bond.i == site0 || bond.j == site0) {
				if(direction0 == 1) {
				//	p0 = p0 ? p0-1 : opsize-1;
				}
				return 4*p + 2*((1-direction0)/2) + (bond.j==site0);
			}
		}
	}

	return std::nullopt;
}
double frust::worm_traverse_measure() {
	const auto matelem_func = [](double direction, const site_basis::state &sold, const site_basis::state &snew) -> double {
		if(sold.j == snew.j && sold.m == snew.m) {
			return sold.jdim != snew.jdim;
			return direction*(sold.jdim-snew.jdim);
		}
		return 0;
	};
	
	double sign = measure_sign();
	
	if(noper_ == 0) {
		measure.add("TauZProb1", 0);
		measure.add("TauZProb2", sign);
		return 0;
	}

	uint32_t p0 = operators_.size()*random01();
	int site0 = lat_.sites.size()*random01();
	int direction0 = 1-2*(random01()>0.5);

	auto v0opt = find_worm_measure_start(site0, p0, direction0);
	if(!v0opt) {
		measure.add("TauZProb1", 0);
		measure.add("TauZProb2", sign);
		return 0;
	}
	uint32_t v0 = *v0opt;

	auto op0 = operators_[v0/4];
	const auto &basis0 = lat_.get_uc_site(site0).basis;
	int wormfunc0 = random01()*basis0.worms.size();

	int leg0 = v0%4;
	const auto &ls0 = lat_.get_vertex_data(op0.bond()).get_legstate(op0.vertex());
	const auto &state_old0 = basis0.states[ls0[leg0]];
	const auto &state_new0 = basis0.states[basis0.worms[wormfunc0].action[ls0[leg0]]];
	double matelem0 = matelem_func(-direction0,state_old0,state_new0);

	uint32_t v = v0;
	int wormfunc = wormfunc0;

	double mean{};

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

		assert(vertices_[vstep] != -1);
		bool up = leg_out > 1;

		mean = 0;
		if(matelem0 && site0 == 0 && site_idx == 1) {
			if((up && v/4 <= p0 && p0 < vnext/4) ||
			   (vnext/4 >= v/4 && !up && (vnext/4 <= p0 || p0 < v/4)) ||
			   (!up && vnext/4 <= p0 && p0 < v/4) ||
			   (vnext/4 <= v/4 && up && (v/4 <= p0 || p0 < vnext/4))) {
				const auto &state_old = site_out.basis.states[vd.get_legstate(old_op.vertex())[leg_out]];
				const auto &state_new = site_out.basis.states[vd.get_legstate(op.vertex())[leg_out]];

				double matelem = matelem_func((1-2*up), state_old,state_new);

				if(matelem != 0) {
					print_operators();
					std::cout << fmt::format("v0: {}, 4p0: {}, v: {}, leg_out: {}, vnext: {}, worm: {}->{} sign: {}\n", v0, 4*p0, v, leg_out, vnext, state_old.name, state_new.name, sign);
				}

				mean = sign*matelem*matelem0;
			}
		}
		if(vnext != v0 || wormfunc != wormfunc0) {
			measure.add("TauZProb1", mean);
			measure.add("TauZProb2", 0);
		}

		v = vnext;
	} while(v != v0 || wormfunc != wormfunc0);
	
	measure.add("TauZProb1", 0);
	measure.add("TauZProb2", sign);


	return 0;
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
		
		wormlength /= floor(nworm_);

		if(!is_thermalized()) {
			avgwormlen_ += 0.01*(wormlength-avgwormlen_);
			double target_worms = param.get<double>("wormlength_fraction", 2.0)*noper_/avgwormlen_;
			nworm_ += 0.01*(target_worms-nworm_)+tanh(target_worms-nworm_);
			nworm_ = std::clamp(nworm_, 1., 1.+noper_/2.);
		}
	} else {
		double mean = 0;

		for(int i = 0; i < nworm_; i++) {
			worm_traverse_measure();
		}	
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
		operators_.resize(operators_.size() * 1.5 + 10, opercode::make_identity());
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
	    j_est{lat_, sign},
	    mag_est<1, 1>{lat_, T_, sign},
	    mag_est<1, -1>{lat_, T_, sign},
	    mag_est<-1, 1>{lat_, T_, sign},
	    mag_est<-1, -1>{lat_, T_, sign},
	};

	std::array<bool, 5> flags = {
	    settings_.measure_j,     settings_.measure_mag,     settings_.measure_sxmag,
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

void frust::measure_chirality(double sign) {
	assert(false);
	int opsize = operators_.size();
	int ucsize = lat_.uc.sites.size();

	const auto tauz = [&](opercode op) -> double {
		if(op.identity()) {
			return 0.0;
		}

		const auto &b = lat_.bonds[op.bond()];
		const auto &bi = lat_.get_uc_site(b.i).basis.states;
		const auto &bj = lat_.get_uc_site(b.j).basis.states;
		const auto &vd = lat_.get_vertex_data(op.bond());
		const auto &ls = vd.get_legstate(op.vertex());



		if(bi[ls[0]].m == bi[ls[2]].m && bj[ls[1]].m == bj[ls[3]].m
		   && bi[ls[0]].j == bi[ls[2]].j && bj[ls[1]].j == bj[ls[3]].j
		   && ((bi[ls[0]].jdim != bi[ls[2]].jdim && bj[ls[1]].jdim == bj[ls[3]].jdim)
		   || (bi[ls[0]].jdim == bi[ls[2]].jdim && bj[ls[1]].jdim != bj[ls[3]].jdim))) {
			double weight = vd.get_weight(op.vertex())*vd.get_sign(op.vertex());
			double signi = 1-2*((b.i/(lat_.Lx*ucsize))&1);
			double signj = 1-2*((b.j/(lat_.Lx*ucsize))&1);
			return (signi*(bi[ls[0]].jdim - bi[ls[2]].jdim) + signj*(bj[ls[1]].jdim - bj[ls[3]].jdim))/weight;
		}

		return 0.0;
	};


	int lastp{};
	opercode lastop;	
	double lastfac{};
	for(; lastp < opsize; lastp++) {
		auto op = operators_[lastp];
		double matelem = tauz(op);
		if(!op.identity()) {
			const auto &vd = lat_.get_vertex_data(op.bond());
			lastfac = matelem/vd.get_weight(op.vertex())*vd.get_sign(op.vertex());
			lastop = op;
			break;
		}
	}

	double mean = 0;

	int p0 = lastp;
	for(int ip = 0; ip < opsize; ip++) {
		int p = (p0+1+ip) % opsize;
		auto op = operators_[p];
		if(op.identity()) {
			continue;
		}
		double matelem = tauz(op); 

		if(lastfac != 0 && matelem != 0) {
			mean += matelem*lastfac;
		}

		lastfac = matelem;

		lastop = op;
		lastp = p;
	}

	measure.add("SignTauZ", -sign*T_*T_*(noper_-1)*mean/lat_.sites.size()/lat_.sites.size());
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

void frust::register_evalables(loadl::evaluator &eval, const loadl::parser &p) {
	auto unsign = [](const std::vector<std::vector<double>> &obs) {
		return std::vector<double>{obs[0][0]/obs[1][0]};
	};

	double T = p.get<double>("T");
	measurement_settings settings{p};
	lattice lat = lattice_from_param(p, false);

	if(settings.measure_j) {
		j_est{lat,0}.register_evalables(eval);
	}
	
	if(settings.measure_mag) {
		mag_est<1,1>{lat, T, 0}.register_evalables(eval);
	}

	if(settings.measure_sxmag) {
		mag_est<-1,1>{lat, T, 0}.register_evalables(eval);
	}
	
	if(settings.measure_sxsymag) {
		mag_est<1,-1>{lat, T, 0}.register_evalables(eval);
	}

	if(settings.measure_sxsymag) {
		mag_est<-1,-1>{lat, T, 0}.register_evalables(eval);
	}

	if(settings.measure_chirality) {
		eval.evaluate("TauZ", {"TauZProb1", "TauZProb2"}, unsign);
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
