#pragma once
#include <string>

#include "Aster/physics/body.h"
#include "Aster/simulations/sim_obj.h"
#include "Aster/simulations/3d_sim_obj.h"


typedef std::string str; 

namespace Aster{

    size_t body_size = sizeof(Body);

    Simulation* bake(sim_meta _s);
    Simulation3d* bake3d(sim3d_meta _s);
}
