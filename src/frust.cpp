#include "frust.h"
#include <algorithm>
#include "vertex_data.h"
#include "latticedef.h"
#include "mag_est.h"

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
		const auto &b = lat_.bonds[op.bond()];
		std::cout << fmt::format("{}-{}: {}\n", b.i, b.j, op.name(lat_.sites[b.i], lat_.sites[b.j]));
	}
	int idx{};
	for(auto s : spin_) {
		std::cout << s.name(lat_.sites[idx]) << ", ";
	}
	std::cout << "\n";
}

void frust::init() {
	spin_.resize(lat_.sites.size());
	std::transform(lat_.sites.begin(), lat_.sites.end(), spin_.begin(), [&](const auto& st) {
		return jm{static_cast<uint8_t>(jm::local_basis_size(st)*random01())};
	});

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
	int wormfunc0 = random01()*vertex_data::wormfunc_count;

	if(lat_.vertex_transition(op0, v0%4, wormfunc0).invalid()) {
		return 0;
	}

	uint32_t v = v0;
	int wormfunc = wormfunc0;
	static int idx = 0;
	idx++;

	do {
		
		wormlength++;

		auto &op = operators_[v/4];
		auto oldop = op;
		int leg_in = v%4;
		const auto &trans = lat_.vertex_transition(op, leg_in, wormfunc);


		auto [leg_out, wormfunc_out, new_vertex_idx] = trans.scatter(random01());
		op = lat_.vertex_idx_opercode(op.bond(), new_vertex_idx);
		if(wormlength >= noper_) {
			const auto &bond = lat_.bonds[op.bond()];
			const auto &si = lat_.sites[bond.i];
			const auto &sj = lat_.sites[bond.j];
			//std::cout << fmt::format("{}-{}-{} {}: {} {} -> {} {}: {}\n", v/4, idx, oldop.bond(), oldop.name(si, sj), leg_in, vertex_data::wormfuncs[wormfunc].name, leg_out, vertex_data::wormfuncs[wormfunc_out].name, op.name(si,sj));
		}

		uint32_t vstep = 4*(v/4)+leg_out;
		if(vstep == v0 && wormfunc_out == vertex_data::wormfuncs[wormfunc0].inverse_idx) {
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
	for(int i = 0; i < nworm_; i++) {
		wormlength += worm_traverse();
	}
	if(is_thermalized()) {
		measure.add("wormlength_fraction", wormlength/noper_);
	}
	
	wormlength /= nworm_;

	if(!is_thermalized()) {
		avgwormlen_ += 0.01*(wormlength-avgwormlen_);
		double target_worms = param.get<double>("wormlength_fraction", 1.8)*noper_/avgwormlen_;
		nworm_ += 0.01*(target_worms-nworm_)+tanh(target_worms-nworm_);
		nworm_ = std::clamp(nworm_, 1., 1.+noper_/2.);
	}

	for(size_t i = 0; i < spin_.size(); i++) {
		if(v_first_[i] < 0) {
			spin_[i] = jm{static_cast<uint8_t>(jm::local_basis_size(lat_.sites[i]) * random01())};
		} else {
			uint32_t v = v_first_[i];
			auto op = operators_[v/4];
			int leg = v%4;
			spin_[i] = op.leg_state(leg);
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

	for(auto &op : operators_) {
		double p_make_bond = lat_.bonds.size()*1./T_/(operators_.size()-noper_);
		double p_remove_bond = (operators_.size()-noper_+1)*T_/lat_.bonds.size();

		if(op.identity()) {
			int bond = random01() * lat_.bonds.size();
			const auto &b = lat_.bonds[bond];
			jm si = spin_[b.i];
			jm sj = spin_[b.j];

			auto newop = opercode::make_vertex(bond, si, sj, si, sj);
			
			double weight = lat_.vertex_weight(newop);

			if(random01() < p_make_bond * weight) {
				op = newop;
				noper_++;
			}
		} else {
			int bond = op.bond();
			const auto &b = lat_.bonds[bond];
			if(op.diagonal()) {
				double weight = lat_.vertex_weight(op);
				if(random01()*weight < p_remove_bond) {
					op = opercode::make_identity();
					noper_--;
				}
			} else {
				spin_[b.i].apply(op.action(0));
				spin_[b.j].apply(op.action(1));
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
			spin_[bond.i].apply(op.action(0));
			spin_[bond.j].apply(op.action(1));
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
			sign += lat_.vertex_sign(op) < 0;
		}
	}
	return (sign&1) ? -1 : 1;
}

void frust::do_measurement() {
	double sign = measure_sign();

	opstring_measurement(mag_est<true>{lat_,T_,sign},mag_est<false>{lat_, T_, sign});
	measure.add("Sign", sign);
	measure.add("nOper", static_cast<double>(noper_));

	measure.add("SignEnergy", sign*(-static_cast<double>(noper_) * T_ - lat_.energy_offset) / lat_.spinhalf_count);
}

void frust::checkpoint_write(const loadl::iodump::group &out) {
	std::vector<unsigned int> saveops(operators_.size());
	std::transform(operators_.begin(), operators_.end(), saveops.begin(), [](opercode op) {
		return op.code();
	});
	
	std::vector<unsigned int> savespins(spin_.size());
	std::transform(spin_.begin(), spin_.end(), savespins.begin(), [](jm s) {
		return s.code();
	});
	
	out.write("noper", noper_);
	out.write("avgwormlen_", avgwormlen_);
	out.write("avgwormlen", avgwormlen_);
	out.write("nworm", nworm_);
	out.write("operators", saveops);
	out.write("spin", savespins);
}

void frust::checkpoint_read(const loadl::iodump::group &in) {
	std::vector<unsigned int> saveops;
	std::vector<unsigned int> savespins;

	in.read("noper", noper_);
	in.read("avgwormlen_", avgwormlen_);
	in.read("nworm", nworm_);
	in.read("operators", saveops);

	operators_.resize(saveops.size());
	std::transform(saveops.begin(), saveops.end(), operators_.begin(), [](uint32_t c) { return opercode{c}; });
	

	in.read("spin", savespins);
	spin_.resize(savespins.size());
	std::transform(savespins.begin(), savespins.end(), spin_.begin(), [](uint8_t c) { return jm{c}; });
}

void frust::register_evalables(loadl::evaluator &eval) {
	auto unsign = [](const std::vector<std::vector<double>> &obs) {
		return std::vector<double>{obs[0][0]/obs[1][0]};
	};
	
	mag_est<true>{lat_, T_, 0}.register_evalables(eval);
	mag_est<false>{lat_, T_, 0}.register_evalables(eval);

	eval.evaluate("Energy", {"SignEnergy", "Sign"}, unsign);
}
