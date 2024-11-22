#pragma once

#include "cluster_magnet.h"
#include <loadleveller/loadleveller.h>

std::unique_ptr<cluster_magnet> cluster_magnet_from_param(const loadl::parser &p);
