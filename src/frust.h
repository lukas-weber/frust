#pragma once

#include <cstdint>
#include <loadleveller/loadleveller.h>
#include <optional>
#include <vector>

#include "sse_data.h"
#include "models/model.h"
#include "measurement_settings.h"
#include "opercode.h"

class frust : public loadl::mc {
public:
	void init() override;
	void do_update() override;

	void do_measurement() override;

	void checkpoint_write(const loadl::iodump::group &out) override;
	void checkpoint_read(const loadl::iodump::group &in) override;

	static void register_evalables(loadl::evaluator &eval, const loadl::parser &p);
	double pt_weight_ratio(const std::string &param_name, double new_param) override;
	void pt_update_param(const std::string &param_name, double new_param) override;

	frust(const loadl::parser &p);
private:
	static constexpr int dump_version_ = 2;
	double T_{};

	int maxwormlen_{0};

	double avgwormlen_{1};
	double nworm_{5};
	int64_t noper_{};

	std::vector<opercode> operators_;
	std::vector<state_idx> spin_;

	const std::unique_ptr<model> model_;
	const sse_data data_;

	measurement_settings settings_;

	void diagonal_update();

	void make_vertex_list();

	bool worm_too_long(int wormlen) const;
	int worm_traverse();
	std::optional<int32_t> find_worm_measure_start(int site0, int32_t &p0, int direction0) const;
	int worm_traverse_measure(double &sign, std::vector<double> &corr);
	bool worm_update();

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
