#pragma once

#include <loadleveller/loadleveller.h>
#include "cluster_magnet.h"

std::unique_ptr<cluster_magnet> cluster_magnet_from_param(const loadl::parser &p);

