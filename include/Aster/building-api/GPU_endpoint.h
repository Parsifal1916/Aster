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
    cl_mem create_vbuffer(std::vector<REAL> base);
    cl_mem create_vbuffer(std::vector<vec3> base);
    cl_mem create_vbuffer(std::vector<cl_ulong> base);
    cl_mem create_vbuffer(std::vector<cl_int> base);
    cl_mem create_vbuffer(size_t bytes);
    void reset_GPU();
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
    cl_kernel compile_kernel(std::string* name, std::string* source, REAL softening, bool del = true);
    
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
}
}