#pragma once

#include "cavity_magnet.h"
#include <loadleveller/loadleveller.h>

std::unique_ptr<cavity_magnet> cavity_magnet_from_param(const loadl::parser &p);
