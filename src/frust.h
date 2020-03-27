#pragma once

#include <loadleveller/loadleveller.h>
#include <vector>
#include <cstdint>

#include "opercode.h"
#include "lattice.h"

class frust : public loadl::mc {
private:
	double T_{};

	uint64_t noper_{};
	std::vector<opercode> operators_;
	std::vector<jm> spin_;

	std::unique_ptr<lattice> lat_;

	void diagonal_update();

	void make_vertex_list();
	void worm(int v0);
	void worm_update();

	void print_vertices();
	void print_operators();

	// temporary vertex list stuff
	std::vector<int32_t> vertices_;
	std::vector<int> v_first_;
	std::vector<int> v_last_;
public:
	void init() override;
	void do_update() override;

	void do_measurement() override;
	void measure_bond_correlations();

	template<class... Estimators>
	void opstring_measurement(Estimators... est);

	void checkpoint_write(const loadl::iodump::group &out) override;
	void checkpoint_read(const loadl::iodump::group &in) override;

	void register_evalables(loadl::evaluator &eval) override;

	frust(const loadl::parser &p);
};

