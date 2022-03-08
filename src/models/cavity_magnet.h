#pragma once

#include "model.h"
#include "util/lattice.h"

class cavity_magnet : public model {
public:
	struct bond {
		double J{1};
	};
	struct mode {
		double omega{};
		double coupling{};
		int max_bosons{};
	};
		
	lattice lat;

	cavity_magnet(const lattice &lat, const std::vector<mode>& modes, const std::vector<bond>& bonds);

	sse_data generate_sse_data() const override;
	void to_json(nlohmann::json& out) const override;
	int normalization_site_count() const override {
		return lat.site_count();
	}
private:
	std::vector<mode> modes_;
	std::vector<bond> bonds_;
};
