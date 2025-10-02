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

#include "Aster/simulations/barnes-hut.h"

namespace Aster{   
namespace Barnes{



/*                   //===---------------------------------------------------------===//
.                    // Barnes-Hut GPU definition                                     //
.                    //===---------------------------------------------------------===*/


template <typename T>
class BHG final : public Barnes_Hut<T> {
    public: 

    BHG();
    void step() override;
    Simulation<T>* load() override;

    ~BHG();

    private:
    cl_kernel tree_builder; 
    cl_kernel force_calculator;     

    void compose_force_kernel();
    void build_tree();
    void upload_force_calc(int num_leaves, int tree_size, int num_bodies);
    void sort_mortons(std::vector<std::pair<uint32_t, uint32_t>>& mortons);
    void load_tree_kernel(int n, int num_leaves, void* mortons);
    
};


}
}

#include "Aster/impl/barnes_hut_GPU.tpp"