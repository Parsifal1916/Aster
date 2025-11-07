#define CL_TARGET_OPENCL_VERSION 300
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstring>
#include <iomanip>

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif

#include "Aster/building-api/logging.h"
#include "Aster/simulations/basic.h"
#include "Aster/simulations/sim_obj.h"

#include "Aster/kernels/eulerian_update.cl.h"
#include "Aster/kernels/generalized_intgrt.cl.h"
#include "Aster/kernels/generalized_intgrt3d.cl.h"
#include "Aster/kernels/newton_update.cl.h"
#include "Aster/kernels/newton_update3d.cl.h"
#include "Aster/kernels/PNs.cl.h"
#include "Aster/kernels/PNs3d.cl.h"

#include "Aster/building-api/GPU_endpoint.h"

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

    cl_ulong global_mem_size;
    clGetDeviceInfo(device, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(global_mem_size), &global_mem_size, NULL);

    log_info(std::string("Available memory on selected device: " + std::to_string(global_mem_size/(1024*1024)) + std::string("MB")));

    
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

inline std::vector<std::pair<cl_program, cl_kernel>> programs;

/**
* @brief compiles a kernel
* @param name: pointer to the name of the function
* @param source: source code of the kernel
* @param k: kernel object onto which to write the kernel 
*/
cl_kernel compile_kernel(std::string* name, std::string* source, REAL softening, bool del) {
    cl_int programResult;
    const char* c = source->c_str();
    size_t l = source->size();

    using clock = std::chrono::high_resolution_clock;
    using ms = std::chrono::duration<double, std::milli>;

    auto start = clock::now();

    // Create the OpenCL program
    cl_program program = clCreateProgramWithSource(context, 1, &c, &l, &programResult);
    if (critical_if(programResult != CL_SUCCESS,
                   "Could not create program: " + std::to_string(programResult)))
        exit(1);

    // Compiler optimization flags
    std::ostringstream opts;
    opts << "-DSOFTENING=" << std::scientific << softening;

    std::string build_options = opts.str();

    cl_device_fp_config fp_config;
    if (clGetDeviceInfo(device, CL_DEVICE_DOUBLE_FP_CONFIG, sizeof(fp_config), &fp_config, nullptr) == CL_SUCCESS) {
        if (fp_config > 0)
            build_options += " -cl-kernel-arg-info";
    }

    cl_int build_result = clBuildProgram(program, 1, &device, build_options.c_str(), nullptr, nullptr);

    size_t log_size = 0;
    clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, nullptr, &log_size);
    if (log_size > 2) {
        std::vector<char> log(log_size);
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, log_size, log.data(), nullptr);
        log_info("Build log for kernel \"" + *name + "\":\n" + std::string(log.begin(), log.end()));
    }

    if (critical_if(build_result != CL_SUCCESS, "Could not build kernel (error " + std::to_string(static_cast<int>(build_result)) + ")")){
        
        std::cout << build_result << " "  << *name;exit(3);
    }

    cl_uint num_ks = 0;
    cl_int enumResult = clCreateKernelsInProgram(program, 0, nullptr, &num_ks);
    if (critical_if(num_ks == 0, "No kernels found in program"))
        exit(-1);

    log_info("Creating kernel \"" + *name + "\"... ", "");

    cl_int kernel_res;
    cl_kernel k = clCreateKernel(program, name->c_str(), &kernel_res);
    if (critical_if(kernel_res != CL_SUCCESS,
                   "Could not load kernel '" + *name + "', error code: " + std::to_string(kernel_res)))
        exit(4);
        
    auto end = clock::now();
    
    if (error_level == LOW_t)
    std::cout << std::setfill(' ') << std::setw(37 - name -> size() ) << "done (" << std::fixed
                << std::setprecision(3)
                << std::setw(time_characters) << ms(end - start).count() / 1000 << " s)\n";
    
    if (del){
        programs.emplace_back(program, k);
    } else {
        clRetainProgram(program);
        clRetainKernel(k);
    }
    return k;
}

inline std::vector<cl_mem> buffers; 

cl_mem create_vbuffer(std::vector<REAL> base){
    size_t bytes = sizeof(REAL) * base.size();
    cl_int err;
    cl_mem ret = clCreateBuffer(context, CL_MEM_READ_WRITE, bytes, NULL, &err);
    Check(err);
    Check(clEnqueueWriteBuffer(queue, ret,  CL_TRUE, 0, bytes, base.data(), 0, NULL, NULL));
    buffers.push_back(ret);
    return ret;
}

cl_mem create_vbuffer(std::vector<vec3> base){
    size_t bytes = sizeof(vec3) * base.size();
    cl_int err;
    cl_mem ret = clCreateBuffer(context, CL_MEM_READ_WRITE, bytes, NULL, &err);
    Check(err);
    Check(clEnqueueWriteBuffer(queue, ret,  CL_TRUE, 0, bytes, base.data(), 0, NULL, NULL));
    buffers.push_back(ret);
    return ret;
}

cl_mem create_vbuffer(std::vector<cl_ulong> base){
    size_t bytes = sizeof(cl_ulong) * base.size();
    cl_int err;
    cl_mem ret = clCreateBuffer(context, CL_MEM_READ_WRITE, bytes, NULL, &err);
    Check(err);
    Check(clEnqueueWriteBuffer(queue, ret,  CL_TRUE, 0, bytes, base.data(), 0, NULL, NULL));
    buffers.push_back(ret);
    return ret;
}

cl_mem create_vbuffer(std::vector<cl_int> base){
    size_t bytes = sizeof(cl_int) * base.size();
    cl_int err;
    cl_mem ret = clCreateBuffer(context, CL_MEM_READ_WRITE, bytes, NULL, &err);
    Check(err);
    Check(clEnqueueWriteBuffer(queue, ret,  CL_TRUE, 0, bytes, base.data(), 0, NULL, NULL));
    buffers.push_back(ret);
    return ret;
}


cl_mem create_vbuffer(size_t bytes){
    cl_int err;
    cl_mem ret = clCreateBuffer(context, CL_MEM_READ_WRITE, bytes, NULL, &err);
    Check(err);
    buffers.push_back(ret);
    return ret;
}

void reset_GPU(){
    clFinish(queue);
    for (auto prog : programs){
        clReleaseProgram(prog.first);
        clReleaseKernel(prog.second);
    }
    for (auto buff : buffers){
        clReleaseMemObject(buff);
    }

    buffers.clear();
    programs.clear();
}


//===---------------------------------------------------------===//
// Bodies Update (UB) compiling                                  //
//===---------------------------------------------------------===//



struct uf_functor { 
    cl_kernel k;
    void operator()(Simulation* _s){
        upload_force_kernel(k, _s);
    }
};


struct ub_functor {
    cl_kernel k;
    int order;
    void operator()(Simulation* _s){
        static const REAL* sequence = saba_coeffs[order];
        // fstep
        for (int i = 0; i < saba_coeff_lng[order]; i+=2){
            upload_update_kernel(k, _s, sequence[i+1], sequence[i]);
            _s -> solver -> compute_forces();
        }
    }
};

func_ptr compile_ub_saba(int ord) {
    log_info("compiling GPU scripts from source...");
    
    static std::string* update_kernel = &general_saba3d;
    static std::string kernel_name = "saba";

    cl_kernel k;

    if (critical_if(ord < 0 || ord > 9, 
    "the GPU-accelerated function \"" + kernel_name + "\" only supports orders from 0 - 9"))
        exit(1);
    

    log_info("done compiling updating scripts");
    return ub_functor{compile_kernel(&kernel_name, update_kernel, 0), ord};
}


//===---------------------------------------------------------===//
// Force Updaters (UF) compiling                                 //
//===---------------------------------------------------------===//


func_ptr compile_uf(force_type t, REAL softening) {
    log_info("compiling GPU scripts from source...");
    size_t index = static_cast<size_t>(t);

    static std::string force_kernels2d[] = {newton_cl, cl_pn1, cl_pn2, cl_pn25};
    static std::string force_kernels3d[] = {newton_cl3d, cl3d_pn1, cl3d_pn2, cl3d_pn25};

    std::string& force_kernel = force_kernels3d[index];

    std::string kernel_names[] = {"cl3d_newton", "pn1", "pn2", "pn25"};

    if (critical_if(index > 3 || index < 0, 
    "could not find a suitable GPU-accelerated function for kernel " + kernel_names[index]))
        exit(1);

    log_info("done compiling force calculation scripts");
    return uf_functor{compile_kernel(&kernel_names[index], &force_kernel, softening)};
}

}
}