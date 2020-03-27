#include "frust.h"

frust::frust(const loadl::parser &p)
	: loadl::mc(p) {
	T_ = param.get<double>("T");
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
		const auto &b = lat_->bonds[op.bond()];
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
	if(noper_ > operators_.size() * 0.5) {
		operators_.resize(operators_.size() * 1.5 + 10, opercode::make_identity());
	}

	for(auto &op : operators_) {
		double p_make_bond = lat_->bonds.size()*1./T_/(operators_.size()-noper_);
		double p_remove_bond = (operators_.size()-noper_+1)*T_/lat_->bonds.size();
		int bond = random01() * lat_->bonds.size();
		const auto &b = lat_->bonds[bond];
		if(op.identity()) {
			jm si = spin_[b.i];
			jm sj = spin_[b.j];

			auto newop = opercode::make_vertex(bond, si, sj, si, sj);
			
			double weight = lat_->vertex_weight(newop);
			if(newop.identity()) {
				continue;
			}

			if(random01() < p_make_bond * weight) {
				op = newop;
				noper_++;
			}
		} else {
			if(op.diagonal()) {
				double weight = lat_->vertex_weight(op);
				if(random01()*weight < p_remove_bond) {
					op = opercode::make_identity();
				}
			} else {
				spin_[b.i].apply(op.action(0));
				spin_[b.i].apply(op.action(1));
			}
		}
	}
}

void frust::do_update() {
	diagonal_update();
	make_vertex_list();
	worm_update();
}

void frust::do_measurement() {
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
	out.write("operators", saveops);
	out.write("spin", savespins);
}

void frust::checkpoint_read(const loadl::iodump::group &in) {
	std::vector<unsigned int> saveops;
	std::vector<unsigned int> savespins;

	in.read("noper", noper_);
	in.read("operators", saveops);

	operators_.resize(saveops.size());
	std::copy(saveops.begin(), saveops.end(), operators_.begin());
	

	in.read("spin", savespins);
	spin_.resize(savespins.size());
	std::copy(savespins.begin(), savespins.end(), spin_.begin());
}

void frust::register_evalables(loadl::evaluator &eval) {
	
}
