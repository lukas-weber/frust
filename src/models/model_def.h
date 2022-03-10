#pragma once

#include "model.h"
#include <loadleveller/loadleveller.h>

std::unique_ptr<model> model_from_param(const loadl::parser &p);
