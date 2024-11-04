#include "frust.h"
#include "models/cavity_magnet/cavity_magnet.h"
#include "models/cavity_magnet/opstring_estimators.h"
#include "models/cluster_magnet/opstring_estimators.h"
#include "models/model_def.h"
#include "vertex_data.h"
#include <algorithm>

#include "template_selector.h"

frust::frust(const loadl::parser &p)
    : loadl::mc(p), model_{model_from_param(p)}, data_{model_->generate_sse_data()} {
	T_ = param.get<double>("T");
	v_first_.resize(data_.site_count);
	v_last_.resize(data_.site_count);

	nworm_ = param.get<double>("init_num_worms", 5);

	maxwormlen_ = param.get<int>("maxwormlen", 0);
	if(maxwormlen_ != 0) {
		std::cout << "Warning: maxwormlen set: this feature is not tested!\n";
	}

	if(param.get<bool>("print_vertex_data", false)) {
		data_.print();
	}
}

void frust::print_operators() {
	int p = -1;
	for(auto op : operators_) {
		p++;
		if(op.identity()) {
			continue;
		}
		const auto &b = data_.get_bond(op.bond());
		const auto &ls = data_.get_vertex_data(op.bond()).get_legstate(op.vertex());
		std::cout << fmt::format("{} {}: {}{}->{}{}\n", data_.nlegs * p,
		                         fmt::join(b, b + data_.nlegs / 2, "-"), ls[0], ls[1], ls[2],
		                         ls[3]);
	}

	std::cout << fmt::format("{}\n", fmt::join(spin_.begin(), spin_.end(), ", "));
}

void frust::print_vertices() {
	for(int p = 0; p < static_cast<int>(operators_.size()); p++) {
		if(vertices_[data_.nlegs * p] > 0) {
			std::cout << fmt::format("{:4d}| {:4d}\n", data_.nlegs * p,
			                         fmt::join(&vertices_[data_.nlegs * p],
			                                   &vertices_[data_.nlegs * p] + data_.nlegs, ","));
		}
	}
	std::cout << fmt::format("vfirst: {}\n", fmt::join(v_first_.begin(), v_first_.end(), ", "));
}

void frust::init() {
	spin_.resize(data_.site_count);
	for(size_t i = 0; i < spin_.size(); i++) {
		spin_[i] = data_.get_site_data(i).dim * random01();
	}

	// FIXME: make this (not) depend on the energy scale of the model
	operators_.resize(param.get("init_opstring_cutoff", static_cast<int>(data_.site_count * T_)));

	const int warmup = 5;
	for(int i = 0; i < warmup; i++) {
		diagonal_update();
	}
}

bool frust::worm_too_long(int wormlen) const {
	if(maxwormlen_ != 0 && wormlen >= maxwormlen_ * static_cast<int>(operators_.size()))
		return true;
	else
		return false;
}

// trick 17: avoid idivs by pulling this kind of hinting trick. These were significant in
// benchmarks.
inline void hint_leg_count(int nlegs, const std::function<void(int)> &func) {
	if(nlegs == 4) {
		func(4);
	} else if(nlegs == 6) {
		func(6);
	} else {
		assert(false);
		__builtin_unreachable();
	}
}

int frust::worm_traverse() {
	if(noper_ == 0) {
		return 0;
	}

	int wormlength{1};
	hint_leg_count(data_.nlegs, [&](int nlegs) {
		int v0{};
		do {
			v0 = vertices_.size() * random01();
		} while(vertices_[v0] < 0);

		auto op0 = operators_[v0 / nlegs];
		int site0 = data_.get_bond(op0.bond())[v0 % (nlegs / 2)];
		int wormfunc0 = random01() * worm_count(data_.get_site_data(site0).dim);

		int v = v0;
		int wormfunc = wormfunc0;

		do {
			auto &op = operators_[v / nlegs];
			assert(!op.vertex().invalid());
			int leg_in = v % nlegs;
			const auto [leg_out, wormfunc_out, new_vertex] =
			    data_.get_vertex_data(op.bond()).scatter(op.vertex(), leg_in, wormfunc, random01());

			op = opercode{op.bond(), new_vertex};

			int vstep = nlegs * (v / nlegs) + leg_out;
			const auto &site_out =
			    data_.get_site_data(data_.get_bond(op.bond())[leg_out % (nlegs / 2)]);
			if((vstep == v0 && wormfunc_out == worm_inverse(wormfunc0, site_out.dim)) ||
			   worm_too_long(wormlength)) {
				break;
			}
			wormlength++;

			wormfunc = wormfunc_out;
			v = vertices_[vstep];
			assert(vertices_[vstep] != -1);
		} while((v != v0 || wormfunc != wormfunc0) && !(worm_too_long(wormlength)));
	});
	return wormlength;
}

std::optional<int> frust::find_worm_measure_start(int site0, int &p0, int direction0) const {
	int opsize = operators_.size();

	if(v_first_[site0] < 0) {
		return std::nullopt;
	}

	for(int l = 0; l <= opsize; l++) {
		int p = p0 + direction0 * l;
		if(p < 0) {
			p += opsize;
		}
		if(p >= opsize) {
			p -= opsize;
		}

		auto op = operators_[p];
		if(!op.identity() && (l != 0 || direction0 < 0)) {
			const auto &bond = data_.get_bond(op.bond());
			const auto &match_site = std::find(bond, bond + data_.nlegs / 2, site0) - bond;

			if(match_site != data_.nlegs / 2) {
				/*if(direction0 == 1) {
				    p0 = p0 ? p0-1 : opsize-1;
				}*/
				return data_.nlegs * p + data_.nlegs / 2 * (direction0 < 0) + match_site;
			}
		}
	}
	assert(false);

	return std::nullopt;
}

int frust::worm_traverse_measure(double &sign, std::vector<double> &corr) {
	if(model_->type != model::model_type::cluster_magnet) {
		throw std::runtime_error{"not implemented"};
	}
	const auto &cm_model = static_cast<cluster_magnet &>(*model_);
	const auto matelem_idx = [](bool switcheroo, const site_basis &basis, state_idx sbefore,
	                            state_idx safter) -> double {
		const auto &sup = switcheroo ? sbefore : safter;
		const auto &sdown = switcheroo ? safter : sbefore;

		if(basis.j(sdown) == basis.j(sup) && basis.m(sdown) == basis.m(sup)) {
			return basis.jdim(sup);
		}
		return -1;
	};

	// assert(data_.uc.sites.size() == 1); // not implemented

	int wormlength{};

	if(noper_ == 0) {
		return wormlength;
	}

	int p0 = operators_.size() * random01();
	int site0 = data_.site_count * random01();
	int direction0 = 1 - 2 * (random01() > 0.5);

	auto v0opt = find_worm_measure_start(site0, p0, direction0);
	if(!v0opt) {
		return wormlength;
	}
	int v0 = *v0opt;
	int v1 = vertices_[v0];
	auto [uc0, x0, y0] = cm_model.lat.split_idx(site0);

	const auto &basis0 = cm_model.get_site(site0).basis;
	int dim0 = data_.get_site_data(site0).dim;
	int wormfunc0 = random01() * worm_count(dim0);

	int v = v0;
	int wormfunc = wormfunc0;

	do {
		wormlength++;
		auto &op = operators_[v / data_.nlegs];
		assert(!op.vertex().invalid());
		opercode old_op = op;

		int leg_in = v % data_.nlegs;
		const auto &vd = data_.get_vertex_data(op.bond());
		const auto [leg_out, wormfunc_out, new_vertex] =
		    vd.scatter(op.vertex(), leg_in, wormfunc, random01());

		op = opercode{op.bond(), new_vertex};

		int vstep = data_.nlegs * (v / data_.nlegs) + leg_out;
		int site_idx = data_.get_bond(op.bond())[leg_out % (data_.nlegs / 2)];

		const auto &site_out = data_.get_site_data(site_idx);
		const auto &model_site_out = cm_model.get_site(site_idx);

		sign *= vd.get_sign(old_op.vertex()) * vd.get_sign(op.vertex());

		if(vstep == v0 && wormfunc_out == worm_inverse(wormfunc0, site_out.dim)) {
			break;
		}

		wormfunc = wormfunc_out;
		int vnext = vertices_[vstep];

		assert(vnext != -1);

		bool up = leg_out > 1;
		if((up && v / data_.nlegs <= p0 && p0 < vnext / data_.nlegs) ||
		   (!up && vnext / data_.nlegs <= p0 && p0 < v / data_.nlegs) ||
		   (up && vnext / data_.nlegs <= v / data_.nlegs &&
		    (p0 >= v / data_.nlegs || p0 < vnext / data_.nlegs)) ||
		   (!up && v / data_.nlegs <= vnext / data_.nlegs &&
		    (p0 >= vnext / data_.nlegs || p0 < v / data_.nlegs))) {
			int state_after_idx = vd.get_legstate(op.vertex())[leg_out];
			int state_before_idx =
			    worm_action(worm_inverse(wormfunc, site_out.dim), state_after_idx, site_out.dim);

			int matidx = matelem_idx(!up, model_site_out.basis, state_before_idx, state_after_idx);

			if(matidx >= 0 && site_idx != site0) {
				auto op0 = operators_[v0 / data_.nlegs];
				auto op1 = operators_[v1 / data_.nlegs];
				const auto &ls0 = data_.get_vertex_data(op0.bond()).get_legstate(op0.vertex());
				const auto &ls1 = data_.get_vertex_data(op1.bond()).get_legstate(op1.vertex());
				int matidx0 = matelem_idx(direction0 > 0, basis0, ls0[v0 % data_.nlegs],
				                          ls1[v1 % data_.nlegs]);

				//				std::cout << fmt::format("{}, {} {}, {}", matidx0,
				// state_before0.name, state_after0.name, wormfunc0) << "\n";
				if(matidx0 >= 0) {
					auto [uc, x, y] = cm_model.lat.split_idx(site_idx);

					int idx = ((y - y0 + cm_model.lat.Ly) % cm_model.lat.Ly) * cm_model.lat.Lx +
					          (x - x0 + cm_model.lat.Lx) % cm_model.lat.Lx;
					if(cm_model.settings.loopcorr_as_strucfac) {
						idx = 0;
					}
					corr[4 * idx + 2 * matidx0 + matidx] += worm_count(dim0) * sign;
				}
			}
		}

		v = vnext;
	} while(v != v0 || wormfunc != wormfunc0);
	return wormlength;
}

bool frust::worm_update() {
	std::vector<opercode> old_operators = {};
	if(maxwormlen_ != 0) {
		old_operators = operators_;
	}
	if(!is_thermalized() || model_->type != model::model_type::cluster_magnet ||
	   (model_->type == model::model_type::cluster_magnet &&
	    !static_cast<const cluster_magnet &>(*model_).settings.measure_chirality)) {
		double totalwormlength = 1;
		for(int i = 0; i < nworm_; i++) {
			int wormlength{};
			wormlength = worm_traverse();
			totalwormlength += wormlength;
			if(worm_too_long(wormlength)) {
				operators_ = old_operators;
				return 0;
			}
		}

		if(is_thermalized() && noper_ != 0) {
			measure.add("wormlength_fraction", totalwormlength / noper_);
		}

		totalwormlength /= ceil(nworm_);

		if(!is_thermalized()) {
			double num_worms_attenuation_factor = param.get<double>("num_worms_attenuation_factor");
			avgwormlen_ += num_worms_attenuation_factor * (totalwormlength - avgwormlen_);
			double target_worms =
			    param.get<double>("wormlength_fraction", 2.0) * noper_ / avgwormlen_;
			nworm_ += num_worms_attenuation_factor * (target_worms - nworm_ + 100*tanh(target_worms - nworm_));
			if(num_worms_attenuation_factor != 0) {
				nworm_ = std::clamp(nworm_, 1., 1. + noper_ / 2.);
			}
		}
	} else {
		auto cm_model = static_cast<const cluster_magnet &>(*model_);
		int size = cm_model.settings.loopcorr_as_strucfac ? 1 : data_.site_count;
		std::vector<double> corr(4 * size);
		double sign = measure_sign();

		for(int i = 0; i < nworm_; i++) {
			worm_traverse_measure(sign, corr);
		}
		for(auto &c : corr) {
			c *= 1. / ceil(nworm_);
		}

		measure.add("SignJDimOffCorr", corr);
	}

	for(size_t i = 0; i < spin_.size(); i++) {
		if(v_first_[i] < 0) {
			spin_[i] = data_.get_site_data(i).dim * random01();
		} else {
			int v = v_first_[i];
			auto op = operators_[v / data_.nlegs];
			int leg = v % data_.nlegs;
			const auto &leg_state = data_.get_vertex_data(op.bond()).get_legstate(op.vertex());
			spin_[i] = leg_state[leg];
		}
	}

	return 1;
}

void frust::make_vertex_list() {
	vertices_.resize(operators_.size() * data_.nlegs);
	std::fill(vertices_.begin(), vertices_.end(), -1);
	std::fill(v_first_.begin(), v_first_.end(), -1);
	std::fill(v_last_.begin(), v_last_.end(), -1);

	for(size_t p = 0; p < operators_.size(); p++) {
		auto op = operators_[p];
		if(op.identity()) {
			continue;
		}
		int v0 = data_.nlegs * p;
		const auto &b = data_.get_bond(op.bond());
		for(int s = 0; s < data_.nlegs / 2; s++) {
			int v1 = v_last_[b[s]];
			if(v1 != -1) {
				vertices_[v1] = v0 + s;
				vertices_[v0 + s] = v1;
			} else {
				v_first_[b[s]] = v0 + s;
			}
			v_last_[b[s]] = v0 + data_.nlegs / 2 + s;
		}
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
	double p_make_bond_raw = data_.bond_count * 1. / T_;
	double p_remove_bond_raw = T_ / data_.bond_count;

	for(auto &op : operators_) {
		if(op.identity()) {
			uint32_t bond = random01() * data_.bond_count;
			const auto &b = data_.get_bond(bond);
			int state_idx = 0;
			for(int s = 0; s < data_.nlegs / 2; s++) {
				state_idx *= data_.get_site_data(b[s]).dim;
				state_idx += spin_[b[s]];
			}

			vertexcode newvert = data_.get_vertex_data(bond).get_diagonal_vertex(state_idx);
			double weight = data_.get_vertex_data(bond).get_weight(newvert);
			double p_make_bond = p_make_bond_raw / (opersize - noper_);

			if(random01() < p_make_bond * weight) {
				op = opercode{bond, newvert};
				noper_++;
			}
		} else {
			int bond = op.bond();
			const auto &vertdata = data_.get_vertex_data(bond);
			if(op.diagonal()) {
				double weight = vertdata.get_weight(op.vertex());
				double p_remove_bond = (opersize - noper_ + 1) * p_remove_bond_raw;
				if(random01() * weight < p_remove_bond) {
					op = opercode::make_identity();
					noper_--;
				}
			} else {
				const auto &b = data_.get_bond(bond);
				const auto &leg_state = vertdata.get_legstate(op.vertex());
				for(int s = 0; s < vertdata.leg_count / 2; s++) {
					spin_[b[s]] = leg_state[vertdata.leg_count / 2 + s];
				}
			}
		}
	}

	assert(tmpspin == spin_);
}

void frust::do_update() {
	bool worm_success{};
	do {
		diagonal_update();
		make_vertex_list();
		worm_success = worm_update();
	} while(!worm_success);
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
			const auto &b = data_.get_bond(op.bond());
			const auto &vd = data_.get_vertex_data(op.bond());
			const auto &leg_state = vd.get_legstate(op.vertex());
			for(int s = 0; s < vd.leg_count / 2; s++) {
				spin_[b[s]] = leg_state[vd.leg_count / 2 + s];
			}
		}

		if(n < noper_) {
			auto x2 = {(ests.measure(op, spin_, data_), 0)...};
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
			sign += data_.get_vertex_data(op.bond()).get_sign(op.vertex()) < 0;
		}
	}
	return (sign & 1) ? -1 : 1;
}

void frust::do_measurement() {
	double sign = measure_sign();

	if(model_->type == model::model_type::cluster_magnet) {
		const auto &cm_model = static_cast<const cluster_magnet &>(*model_);
		auto [obs, flags] = cluster_magnet_opstring_estimators(cm_model, T_, sign);
		template_select([&](auto... vals) { opstring_measurement(vals...); }, obs, flags);
	}

	if(model_->type == model::model_type::cavity_magnet) {
		const auto &cm_model = static_cast<const cavity_magnet &>(*model_);
		auto [obs, flags] = cavity_magnet_opstring_estimators(cm_model, T_, sign);
		template_select([&](auto... vals) { opstring_measurement(vals...); }, obs, flags);
	}

	measure.add("Sign", sign);
	measure.add("nOper", static_cast<double>(noper_));
	measure.add("SignNOper", sign * static_cast<double>(noper_));
	measure.add("SignNOper2", sign * static_cast<double>(noper_) * static_cast<double>(noper_));

	measure.add("SignEnergy", sign * (-static_cast<double>(noper_) * T_ - data_.energy_offset) /
	                              model_->normalization_site_count());
}

void frust::checkpoint_write(const loadl::iodump::group &out) {
	std::vector<opercode_uint> saveops(operators_.size());
	std::transform(operators_.begin(), operators_.end(), saveops.begin(),
	               [](opercode op) { return op.code(); });

	out.write("noper", noper_);
	out.write("avgwormlen", avgwormlen_);
	out.write("nworm", nworm_);
	out.write("operators", saveops);
	out.write("spin", spin_);
	out.write("version", dump_version_);
}

void frust::checkpoint_read(const loadl::iodump::group &in) {
	std::vector<opercode_uint> saveops;

	in.read("noper", noper_);
	in.read("avgwormlen", avgwormlen_);
	in.read("nworm", nworm_);
	in.read("operators", saveops);

	operators_.resize(saveops.size());
	std::transform(saveops.begin(), saveops.end(), operators_.begin(),
	               [](opercode_uint c) { return opercode{c}; });

	in.read("spin", spin_);

	int version;
	in.read("version", version);
	if(version != dump_version_) {
		throw std::runtime_error("dump file does not fit program version");
	}
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
	std::unique_ptr<model> m = model_from_param(p);

	m->register_evalables(eval, T);

	eval.evaluate("Energy", {"SignEnergy", "Sign"}, unsign);
	eval.evaluate("SpecificHeat", {"SignNOper2", "SignNOper", "Sign"},
	              [&](const std::vector<std::vector<double>> &obs) {
		              double sn2 = obs[0][0];
		              double sn = obs[1][0];
		              double sign = obs[2][0];

		              return std::vector<double>{(sn2 / sign - sn * sn / sign / sign - sn / sign) /
		                                         m->normalization_site_count()};
	              });
}

double frust::pt_weight_ratio(const std::string &param_name, double new_param) {
	// we only implement temperature swaps
	if(param_name != "T") {
		return 0;
	}

	return -noper_ * (log(new_param / T_));
}

void frust::pt_update_param(const std::string &param_name, double new_param) {
	if(param_name == "T") {
		T_ = new_param;
	}
}
