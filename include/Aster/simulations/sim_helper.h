#pragma once
#include <string>

#include "Aster/physics/body.h"
#include "Aster/simulations/sim_obj.h"

typedef std::string str; 

namespace Aster{
    Simulation<vec3>* bake3d(simulation_types s);
    Simulation<vec2>* bake(simulation_types s);
}
