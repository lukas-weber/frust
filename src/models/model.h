#pragma once

#include <vector>
#include <nlohmann/json_fwd.hpp>
#include "../sse_data.h"

class model {
public:
	enum class model_type {
		cluster_magnet,
		cavity_magnet
	};

	model_type type;

	model(model_type type) : type{type} {};	
	virtual ~model() = default;
	virtual sse_data generate_sse_data() const = 0;
	virtual void to_json(nlohmann::json& out) const = 0;

	// The site count used for normalization, e.g. in the energy.
	// Some models have multiple physical degrees of freedom per computational
	// site (e.g. cluster_magnet)
	virtual int normalization_site_count() const = 0;
};
