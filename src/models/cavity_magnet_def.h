#pragma once

#include <loadleveller/loadleveller.h>
#include "cavity_magnet.h"

std::unique_ptr<cavity_magnet> cavity_magnet_from_param(const loadl::parser &p);
