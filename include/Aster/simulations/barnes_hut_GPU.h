#pragma once
#include <vector>
#include <cassert>
#include <thread>

#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#define CL_TARGET_OPENCL_VERSION 300

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif

#include "Aster/simulations/sim_obj.h"
#include "Aster/simulations/barnes-hut.h"

namespace Aster{   
namespace Barnes{



/*                   //===---------------------------------------------------------===//
.                    // Barnes-Hut GPU definition                                     //
.                    //===---------------------------------------------------------===*/



class BHG final : public Barnes_Hut {
    public: 

    BHG(Simulation* _s);
    void load() override;
    void compute_forces() override;

    private:
    cl_kernel tree_builder; 
    cl_kernel force_calculator;     

    void compose_force_kernel();
    void build_tree();
    void upload_force_calc(int num_leaves, int tree_size, int num_bodies);
    void load_tree_kernel(int n, int num_leaves, void* mortons);
    
};


}
}