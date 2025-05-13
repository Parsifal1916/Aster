#pragma once
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#define CL_TARGET_OPENCL_VERSION 300

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

namespace Aster{
namespace GPU{


FORCE_INLINE void check(cl_int a){
    static int b = 0;
    b++;
    if (critical_if(a != CL_SUCCESS, "error in calculating the force (GPU side)")){
        std::cout << b << " " << a << "\n";
        exit(-1);
    }
}
    

template <typename T> 
void upload_update_kernel(cl_kernel& k, Simulation<T>* _s, double c, double d){
    if (!has_initialized) {
        init_opencl();
        has_initialized = true;
    }
    
    cl_int operation_result;
    uint8_t vec_size = sizeof(T)/8;


    if (_s -> bodies.positions.size() == 0) return;
    
    const size_t num_bytes = _s -> bodies.positions.size() * sizeof(double);
    const cl_uint N = _s -> bodies.positions.size();
    const size_t arr_size = num_bytes*vec_size;
    const double dt = _s -> get_dt();

    auto pos_b = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, arr_size, _s -> bodies.positions.data(),&operation_result);
    check(operation_result);
    auto vel_b = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, arr_size, _s -> bodies.velocities.data(),&operation_result);
    check(operation_result);
    auto acc_b = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, arr_size, _s -> bodies.accs.data(),&operation_result);
    check(operation_result);

    size_t LW_size = 64;
    size_t GW_size = ((N + LW_size - 1) / LW_size) * LW_size;
    double c_arg = c;
    double d_arg = d;

    operation_result = clSetKernelArg(k, 0, sizeof(unsigned int), &N);
    check(operation_result);
    operation_result = clSetKernelArg(k, 1, sizeof(double), &dt);
    check(operation_result);
    operation_result = clSetKernelArg(k, 2, sizeof(double), &d_arg);
    check(operation_result);
    operation_result = clSetKernelArg(k, 3, sizeof(double), &c_arg);
    check(operation_result);
    operation_result = clSetKernelArg(k, 4, sizeof(cl_mem), &pos_b);
    check(operation_result);
    operation_result = clSetKernelArg(k, 5, sizeof(cl_mem), &vel_b);
    check(operation_result);
    operation_result = clSetKernelArg(k, 6, sizeof(cl_mem), &acc_b);
    check(operation_result);


    operation_result = clEnqueueNDRangeKernel(queue, k, 1, 0, &GW_size, &LW_size, 0, nullptr, nullptr );
    check(operation_result);

    operation_result = clEnqueueReadBuffer(queue, pos_b, CL_TRUE, 0, arr_size, _s -> bodies.positions.data(), 0, nullptr, nullptr);
    check(operation_result);
    operation_result = clEnqueueReadBuffer(queue, vel_b, CL_TRUE, 0, arr_size, _s -> bodies.velocities.data(), 0, nullptr, nullptr);
    check(operation_result);
    operation_result = clEnqueueReadBuffer(queue, acc_b, CL_TRUE, 0, arr_size, _s -> bodies.accs.data(), 0, nullptr, nullptr);
    check(operation_result);
    
    clFinish( queue );

    clReleaseMemObject(pos_b);
    clReleaseMemObject(vel_b);
    clReleaseMemObject(acc_b);
}


template <typename T> 
void upload_force_kernel(cl_kernel& k, Simulation<T>* _s){
    if (!has_initialized) {
        init_opencl();
        has_initialized = true;
    }
    cl_int operation_result;
    if (_s -> bodies.positions.size() == 0) return;

    const size_t num_bytes = _s -> bodies.positions.size() * sizeof(double);
    const cl_uint N = _s -> bodies.positions.size();
    const double G = _s -> get_G();
    const size_t vec_size = sizeof(T)/sizeof(double)*num_bytes;

    // creates buffer for mass
    cl_mem masses_b = clCreateBuffer(context, CL_MEM_READ_ONLY, num_bytes, nullptr, &operation_result);
    check(operation_result);

    // loads buffer for mass
    operation_result = clEnqueueWriteBuffer(queue, masses_b, CL_TRUE, 0, num_bytes, _s -> bodies.masses.data(), 0, nullptr, nullptr);
    check(operation_result);

    // creates buffer for temperature
    cl_mem temp_b = clCreateBuffer(context, CL_MEM_READ_ONLY, num_bytes, nullptr, &operation_result);
    check(operation_result);
    
    // loads buffer for temperature
    operation_result = clEnqueueWriteBuffer(queue, temp_b, CL_TRUE, 0, num_bytes, _s -> bodies.temps.data(), 0, nullptr, nullptr);
    check(operation_result);

    // creates buffer for position
    cl_mem pos_b = clCreateBuffer(context, CL_MEM_READ_ONLY, vec_size, nullptr, &operation_result);
    check(operation_result);

    // loads buffer for position
    operation_result = clEnqueueWriteBuffer(queue, pos_b, CL_TRUE, 0, vec_size, _s -> bodies.positions.data(), 0, nullptr, nullptr);
    check(operation_result);

    // creates buffer for velocities
    cl_mem vel_b = clCreateBuffer(context, CL_MEM_READ_ONLY, vec_size, nullptr, &operation_result);
    check(operation_result);

    // loads buffer for velocities
    operation_result = clEnqueueWriteBuffer(queue, vel_b, CL_TRUE, 0, vec_size, _s -> bodies.velocities.data(), 0, nullptr, nullptr);
    check(operation_result);

    //creates buffer for acceleration
    cl_mem acc_br = clCreateBuffer( context, CL_MEM_WRITE_ONLY, vec_size, nullptr, &operation_result);
    check(operation_result);

    size_t LW_size = 256;
    size_t GW_size = ((N + LW_size - 1) / LW_size) * LW_size;

    operation_result = clSetKernelArg(k, 0, sizeof(unsigned int), &N);
    check(operation_result);
    operation_result = clSetKernelArg(k, 1, sizeof(double), &G);
    check(operation_result);
    operation_result = clSetKernelArg(k, 2, sizeof(cl_mem), &temp_b);
    check(operation_result);
    operation_result = clSetKernelArg(k, 3, sizeof(cl_mem), &masses_b);
    check(operation_result);
    operation_result = clSetKernelArg(k, 4, sizeof(cl_mem), &pos_b);
    check(operation_result);
    operation_result = clSetKernelArg(k, 5, sizeof(cl_mem), &vel_b);
    check(operation_result);
    operation_result = clSetKernelArg(k, 6, sizeof(cl_mem), &acc_br);
    check(operation_result);

    operation_result = clEnqueueNDRangeKernel(queue, k, 1, 0, &GW_size, &LW_size, 0, nullptr, nullptr );
    check(operation_result);

    operation_result = clEnqueueReadBuffer( queue, acc_br, CL_TRUE, 0, vec_size, _s -> bodies.accs.data(), 0, nullptr, nullptr );
    check(operation_result);

    clFinish(queue);

    clReleaseMemObject(pos_b);
    clReleaseMemObject(acc_br);
    clReleaseMemObject(masses_b);
    clReleaseMemObject(temp_b);
    clReleaseMemObject(vel_b);
}


}
}