#pragma once

#include <loadleveller/loadleveller.h>
#include "model.h"

std::unique_ptr<model> model_from_param(const loadl::parser &p);
