#pragma once
#define CL_TARGET_OPENCL_VERSION 300
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS

#include <fstream>
#include <string>
#include <cstring>

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif

#include "Aster/building-api/logging.h"
#include "Aster/simulations/basic.h"
#include "Aster/impl/config.h"
#include "Aster/simulations/sim_obj.h"

#include "Aster/impl/kernels/eulerian_update.cl.h"
#include "Aster/impl/kernels/generalized_intgrt.cl.h"
#include "Aster/impl/kernels/generalized_intgrt3d.cl.h"
#include "Aster/impl/kernels/newton_update.cl.h"
#include "Aster/impl/kernels/newton_update3d.cl.h"
#include "Aster/impl/kernels/PNs.cl.h"
#include "Aster/impl/kernels/PNs3d.cl.h"

#include "Aster/impl/kernel_uploaders.tpp"

namespace Aster{
namespace GPU{

//===---------------------------------------------------------===//
// Initilizing                                                   //
//===---------------------------------------------------------===//

/**
* @brief sets the best device in terms of compute power
*/
void select_best_device() {
    log_info("choosing best device...");
    // gathers all the platforms
    cl_uint num_platforms;
    cl_platform_id platforms[64];
    cl_uint result = clGetPlatformIDs(64, platforms, &num_platforms);
    if (critical_if(result != CL_SUCCESS || num_platforms == 0, "no OpenCL platforms found, you might want to install the correct drivers or OpenCL2")) exit(-1);


    // pre-loop setup
    cl_device_id best_device = nullptr;
    cl_platform_id best_plat = nullptr;
    cl_uint max_score = 0;

    for (cl_uint i = 0; i < num_platforms; ++i) {//goes through every platform
        cl_uint num_devices;
        cl_device_id devices[64];
        clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, 64, devices, &num_devices);

        for (cl_uint j = 0; j < num_devices; ++j) {//iterates over every device
            // gathers compute units and clock frquency
            cl_uint compute_units;
            cl_uint clock_frequency;
            clGetDeviceInfo(devices[j], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &compute_units, nullptr);
            clGetDeviceInfo(devices[j], CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(cl_uint), &clock_frequency, nullptr);

            // computes score
            cl_uint score = compute_units * clock_frequency;
            if (score <= max_score) continue;

            max_score = score;
            best_device = devices[j];
            best_plat = platforms[i];
        }
    }

    platform = best_plat;
    device = best_device;
}

/**
* @brief initializes opencl and finds the right device
*/
void init_opencl(){
    if (has_initialized) return
    log_info("initializing GPU drivers...");
    select_best_device();

    // sets up the context
	cl_int contextResult;
	context = clCreateContext( nullptr, 1, &device, nullptr, nullptr, &contextResult );
	
    if (critical_if( contextResult != CL_SUCCESS, "could not load an OpenCL context" )) exit(-1);

	cl_int commandQueueResult;
	queue = clCreateCommandQueue( context, device, 0, &commandQueueResult );
	if (critical_if( commandQueueResult != CL_SUCCESS, "could not load an OpenCL command queue")) exit(-1);

    has_initialized = true;
}

/**
* @brief compiles a kernel
* @param name: pointer to the name of the function
* @param source: source code of the kernel
* @param k: kernel object onto which to write the kernel 
*/
void compile_kernel(std::string* name, std::string* source, cl_kernel& k){

    cl_int programResult;
    const char* c = source->c_str();
    size_t l = source->size();

    cl_program program = clCreateProgramWithSource(context, 1, &c, &l, &programResult);

    if (critical_if(programResult != CL_SUCCESS, 
                   "could not create program: " + std::to_string(programResult)))
        exit(1);

    cl_int build_result = clBuildProgram(program, 1, &device, "", nullptr, nullptr);

    size_t log_size = 0;
    clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, nullptr, &log_size);

    if (log_size > 1) {
        std::vector<char> log(log_size);
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, log_size, log.data(), nullptr);
        std::string build_log(log.begin(), log.end());
        log_info("Build log for kernel \"" + *name + "\":\n" + build_log);
    }

    if (critical_if(build_result != CL_SUCCESS, 
                   "could not build kernel"))
        exit(3);

    cl_uint num_ks = 0;
    cl_int enumResult = clCreateKernelsInProgram(program, 0, nullptr, &num_ks);

    if (critical_if(num_ks == 0, "no kernels found in program"))
        exit(-1);

    log_info("creating kernel \"" + *name +"\"");

    cl_int kernel_res;
    k = clCreateKernel(program, name->c_str(), &kernel_res);

    if (critical_if(kernel_res != CL_SUCCESS, 
                   "could not load kernel '" + *name + 
                   "', error code: " + std::to_string(kernel_res))) {
        exit(4);
    }
}

//===---------------------------------------------------------===//
// Bodies Update (UB) compiling                                  //
//===---------------------------------------------------------===//


template <typename T>
struct uf_functor {
    cl_kernel k;
    void operator()(Simulation<T>* _s){
        upload_force_kernel<T>(k, _s);
    }
};

template <typename T>
struct ub_functor {
    cl_kernel k;
    int order;
    void operator()(Simulation<T>* _s){
        static const double* sequence = saba_coeffs[order];

        // fstep
        for (int i = 0; i < saba_coeff_lng[order]; i+=2){
            upload_update_kernel(k, _s, sequence[i], sequence[i+1]);
            _s -> update_forces(_s);
        }
    }
};

template <typename T>
func_ptr<T> compile_ub(update_type t) {
    log_info("compiling GPU scripts from source...");
    
    static std::string* update_kernel = std::is_same<T, vec2>() ? &general_saba : &general_saba3d;
    static std::string kernel_name = "saba";
    int index =  static_cast<int>(t);

    cl_kernel k;

    if (critical_if(index < 0 || index > static_cast<int>(CUSTOM_U), 
    "could not find a suitable GPU-accelerated function for kernel " + kernel_name))
        exit(1);
    
    compile_kernel(&kernel_name, update_kernel, k);

    log_info("done compiling updating scripts");
    return ub_functor<T>{k, index};
}


//===---------------------------------------------------------===//
// Force Updaters (UF) compiling                                 //
//===---------------------------------------------------------===//

template <typename T>
func_ptr<T> compile_uf(force_type t) {
    log_info("compiling GPU scripts from source...");
    size_t index = static_cast<size_t>(t);
    
    std::string* force_kernel = std::is_same<T, vec2>::value 
        ? (t == NEWTON ? &newton_cl : &cl_pns)
        : (t == NEWTON ? &newton_cl3d : &cl_pns3d)
    ;

    std::string kernel_names[] = {"newton", "pn1", "pn2", "pn25"};
    cl_kernel k;

    if (critical_if(index >= static_cast<int>(CUSTOM_FU) || index < 0, 
    "could not find a suitable GPU-accelerated function for kernel " + kernel_names[index]))
        exit(1);
    
    compile_kernel(&kernel_names[index], force_kernel, k);

    log_info("done compiling force calculation scripts");
    return uf_functor<T>{k};
}

}
}