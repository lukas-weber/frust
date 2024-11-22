#pragma once

#include "lattice.h"
#include <loadleveller/loadleveller.h>

lattice lattice_from_param(const loadl::parser &p, bool with_vertex_data);
