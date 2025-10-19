#pragma once
#include <vector>
#include <map>
#include <string>
#define CL_TARGET_OPENCL_VERSION 300
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS


#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif

#include "Aster/physics/body.h"
#include "Aster/simulations/basic.h"

namespace Aster{

namespace GPU{
    //===---------------------------------------------------------===//
    // GPU optimized stuff                                           //
    //===---------------------------------------------------------===//
    
    inline cl_platform_id platform;
    inline cl_device_id device;
    
    inline cl_context context;
    inline cl_command_queue queue;
    
    inline bool has_initialized = false;
    
    
    /**
    * @brief sets the best device in terms of compute power
    */
    void select_best_device();
    
    /**
    * @brief compiles a kernel
    * @param name: pointer to the name of the function
    * @param source: source code of the kernel
    * @param k: kernel object onto which to write the kernel 
    */
    cl_kernel compile_kernel(std::string* name, std::string* source, REAL softening);
    
    /**
    * @brief initializes opencl and finds the right device
    */
    void init_opencl();
    
    /**
    * @brief compiles the force program
    */
    
    func_ptr compile_uf(force_type t, REAL softening);
    
    /**
    * @brief compiles the body update program
    */
    
    func_ptr compile_ub_saba(int ord);
    
     
    void upload_force_kernel(cl_kernel& k, Simulation* _s);
    
     
    void upload_update_kernel(cl_kernel& k, Simulation* _s, REAL c = 1, REAL d = 1);
    
    /**
    * @brief merges two already sorted arrays
    * @param left first array
    * @param L first array size
    * @param right second array
    * @param R second array size
    * @param merged ptr to the merged array
    */
    static void merge_arrays(const uint32_t* left, size_t L,const uint32_t* right, size_t R,uint32_t* merged);
    
    cl_program compile_sorting_kernels();
    
    /**
    * @brief sorts an array using the gpu
    * @param input: ptr to the start of the arrya
    * @param output: array onto which to write
    * @param N: size of the array
    */
    void sort(uint64_t* input, uint64_t* output, size_t N);
    }
    
    
}
