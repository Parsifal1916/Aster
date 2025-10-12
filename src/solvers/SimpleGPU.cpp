#include "Aster/simulations/sim_obj.h"
#include "Aster/building-api/builder.h"
#include "Aster/building-api/GPU_endpoint.h"

namespace Aster{

SimpleGPU::SimpleGPU(Simulation* s){
    this -> _s = s;
}

void SimpleGPU::load(){
    this -> load_force_kernel = GPU::compile_uf(this -> _t);
}

void SimpleGPU::compute_forces(){
    this -> load_force_kernel(this -> _s);
}



}