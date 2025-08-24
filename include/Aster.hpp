// *  █████╗ ███████╗████████╗███████╗██████╗ 
// * ██╔══██╗██╔════╝╚══██╔══╝██╔════╝██╔══██╗
// * ███████║███████╗   ██║   █████╗  ██████╔╝
// * ██╔══██║╚════██║   ██║   ██╔══╝  ██╔══██╗
// * ██║  ██║███████║   ██║   ███████╗██║  ██║
// * ╚═╝  ╚═╝╚══════╝   ╚═╝   ╚══════╝╚═╝  ╚═╝

#define CL_TARGET_OPENCL_VERSION 200

#include "Aster/graphics/3d_graphics.h"
#include "Aster/graphics/2d_graphics.h"
#include "Aster/graphics/text_rendering.h"

#include "Aster/physics/body.h"
#include "Aster/physics/vectors.h"
#include "Aster/physics/tool-chain.h"

#include "Aster/graphs/labeling.h"

#include "Aster/building-api/builder.h"
#include "Aster/building-api/presets.h"
#include "Aster/building-api/logging.h"

#include "Aster/simulations/sim_helper.h"
#include "Aster/simulations/sim_obj.h"
#include "Aster/simulations/barnes-hut.h"
#include "Aster/simulations/barnes_hut_GPU.h"
#include "Aster/simulations/simulation.h"