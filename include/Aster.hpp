// *  █████╗ ███████╗████████╗███████╗██████╗ 
// * ██╔══██╗██╔════╝╚══██╔══╝██╔════╝██╔══██╗
// * ███████║███████╗   ██║   █████╗  ██████╔╝
// * ██╔══██║╚════██║   ██║   ██╔══╝  ██╔══██╗
// * ██║  ██║███████║   ██║   ███████╗██║  ██║
// * ╚═╝  ╚═╝╚══════╝   ╚═╝   ╚══════╝╚═╝  ╚═╝

#include "Aster/renderer_endpoint.hpp"
#include "Aster/physics_endpoint.hpp"
#include "Aster/builder_endpoint.hpp"

#include "Aster/building-api/builder.h"
#include "Aster/building-api/presets.h"

#include "Aster/simulations/sim_helper.h"
#include "Aster/simulations/sim_obj.h"
#include "Aster/simulations/barnes-hut.h"
#include "Aster/simulations/3d_sim_obj.h"
#include "Aster/simulations/barnes-hut3d.h"
#include "Aster/simulations/simulation.h"