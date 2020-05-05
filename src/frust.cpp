#include "frust.h"
#include <algorithm>
#include "vertex_data.h"
#include "latticedef.h"
#include "mag_est.h"
#include "j_est.h"

frust::frust(const loadl::parser &p)
	: loadl::mc(p),
	  lat_{lattice_from_param(p)} {
	T_ = param.get<double>("T");
	v_first_.resize(lat_.sites.size());
	v_last_.resize(lat_.sites.size());

	//lat_.vertex_print();
}

void frust::print_operators() {
	for(auto op : operators_) {
		if(op.identity()) {
			continue;
		}
		//const auto &b = lat_.bonds[op.bond()];
		//std::cout << fmt::format("{}-{}: {}\n", b.i, b.j, op.name(lat_.get_uc_site(b.i).basis, lat_.get_uc_site(b.j).basis));
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
}

int frust::worm_traverse() {
	if(noper_ == 0) { // XXX: remove when loop measurements are added
		return 0;
	}
	
	int wormlength{};

	uint32_t v0{};
	do {
		v0 = vertices_.size()*random01();
	} while(vertices_[v0] < 0);

	auto op0 = operators_[v0/4];
	int wormfunc0 = random01()*site_basis::worm_count;

	if(lat_.get_vertex_data(op0.bond()).get_transition(op0.vertex(), v0%4, wormfunc0).invalid()) {
		return 0;
	}

	uint32_t v = v0;
	int wormfunc = wormfunc0;
	static int idx = 0;
	idx++;

	do {
		wormlength++;

		auto &op = operators_[v/4];
		assert(!op.vertex().invalid());
		//auto oldop = op;
		int leg_in = v%4;
		const auto &trans = lat_.get_vertex_data(op.bond()).get_transition(op.vertex(), leg_in, wormfunc);


		auto [leg_out, wormfunc_out, new_vertex] = trans.scatter(random01());
		op = opercode{op.bond(), new_vertex};/*
		if(wormlength >= noper_) {
			const auto &bond = lat_.bonds[op.bond()];
			const auto &si = lat_.sites[bond.i];
			const auto &sj = lat_.sites[bond.j];
			std::cout << fmt::format("{}-{}-{} {}: {} {} -> {} {}: {}\n", v/4, idx, oldop.bond(), oldop.vertex().vertex_idx(), leg_in, wormfunc, leg_out, wormfunc_out, op.vertex().vertex_idx());
		}*/

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

void frust::worm_update() {
	//std::cout << fmt::format("nworm : {}\n", nworm_);
	double wormlength = 1;
	auto tmpoperators = operators_;
	for(int i = 0; i < nworm_; i++) {
		tmpoperators = operators_;
		int wl = worm_traverse();
		if(wl < 0) {
			operators_ = tmpoperators;
		} else {
			wormlength += wl;
		}
	}
	if(is_thermalized()) {
		measure.add("wormlength_fraction", wormlength/noper_);
	}
	
	wormlength /= nworm_;

	if(!is_thermalized()) {
		avgwormlen_ += 0.01*(wormlength-avgwormlen_);
		double target_worms = param.get<double>("wormlength_fraction", 1.0)*noper_/avgwormlen_;
		nworm_ += 0.01*(target_worms-nworm_)+tanh(target_worms-nworm_);
		nworm_ = std::clamp(nworm_, 1., 1.+noper_/2.);
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

	/*std::vector<std::vector<double>> bond_histograms(lat_.bonds.size());
	for(int i = 0; i < lat_.bonds.size(); i++) {
		bond_histograms[i].resize(lat_.get_vertex_data(i).vertex_count());
	}*/

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
			/*uint32_t v = op.vertex().vertex_idx();
			auto &histogram = bond_histograms[bond];
			histogram[v]++;*/
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
	/*
	if(is_thermalized()) {
		measure.add("bond_hist", bond_histograms[0]);
	}*/
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

	opstring_measurement(
			//   mag_est<true>{lat_,T_,sign},
			     mag_est<false>{lat_, T_, sign}
			//     j_est{lat_, sign}
			     );
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
}

void frust::register_evalables(loadl::evaluator &eval) {
	auto unsign = [](const std::vector<std::vector<double>> &obs) {
		return std::vector<double>{obs[0][0]/obs[1][0]};
	};
	
	//mag_est<true>{lat_, T_, 0}.register_evalables(eval);
	mag_est<false>{lat_, T_, 0}.register_evalables(eval);
	//j_est{lat_,0}.register_evalables(eval);

	eval.evaluate("Energy", {"SignEnergy", "Sign"}, unsign);
	eval.evaluate("SpecificHeat", {"SignNOper2", "SignNOper", "Sign"}, [&](const std::vector<std::vector<double>> &obs) {
		double sn2 = obs[0][0];
		double sn = obs[1][0];
		double sign = obs[2][0];

		return std::vector<double>{(sn2/sign - sn*sn/sign/sign - sn/sign)/lat_.spinhalf_count};
	});
}
