#pragma once

#include <loadleveller/loadleveller.h>
#include <vector>
#include <cstdint>

#include "opercode.h"
#include "lattice.h"

class frust : public loadl::mc {
public:
	void init() override;
	void do_update() override;

	void do_measurement() override;

	void checkpoint_write(const loadl::iodump::group &out) override;
	void checkpoint_read(const loadl::iodump::group &in) override;

	void register_evalables(loadl::evaluator &eval) override;

	frust(const loadl::parser &p);
private:
	double T_{};


	double avgwormlen_{1};
	double nworm_{5};
	int64_t noper_{};
	
	std::vector<opercode> operators_;
	std::vector<jm> spin_;

	lattice lat_;

	void diagonal_update();

	void make_vertex_list();
	int worm_traverse();
	void worm_update();

	template<class... Estimators>
	void opstring_measurement(Estimators... est);
	double measure_sign() const;

	void print_vertices();
	void print_operators();

	// temporary vertex list stuff
	std::vector<int32_t> vertices_;
	std::vector<int> v_first_;
	std::vector<int> v_last_;
};

