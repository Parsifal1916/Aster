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
#include "Aster/building-api/GPU_endpoint.h"
#include "Aster/simulations/basic.h" 
#include "Aster/simulations/sim_obj.h"

namespace Aster{
namespace GPU{

 
void upload_update_kernel(cl_kernel& kernel, Simulation* _s, REAL c, REAL d){
    if (!has_initialized) {
        init_opencl();
        has_initialized = true;
    }

    _s -> solver -> compute_forces();
    
    if (_s -> bodies.positions.size() == 0) return;

    const cl_uint N = _s -> bodies.positions.size();
    const REAL dt = _s -> get_dt();

    size_t LW_size = 64;
    size_t GW_size = ((N + LW_size - 1) / LW_size) * LW_size;
    REAL c_arg = c;
    REAL d_arg = d;

    Check(clSetKernelArg(kernel, 0, sizeof(unsigned int), &N));
    Check(clSetKernelArg(kernel, 1, sizeof(REAL), &dt));
    Check(clSetKernelArg(kernel, 2, sizeof(REAL), &d_arg));
    Check(clSetKernelArg(kernel, 3, sizeof(REAL), &c_arg));
    Check(clSetKernelArg(kernel, 4, sizeof(cl_mem), &_s -> positions_cl));
    Check(clSetKernelArg(kernel, 5, sizeof(cl_mem), &_s -> velocities_cl));
    Check(clSetKernelArg(kernel, 6, sizeof(cl_mem), &_s -> accs_cl));

    Check(clEnqueueNDRangeKernel(queue, kernel, 1, 0, &GW_size, &LW_size, 0, nullptr, nullptr ));
}


 
void upload_force_kernel(cl_kernel& k, Simulation* _s){
    if (!has_initialized) {
        init_opencl();
        has_initialized = true;
    }
    cl_int operation_result;
    if (_s -> bodies.positions.size() == 0) return;

    const size_t num_bytes = _s -> bodies.positions.size() * sizeof(REAL);
    const cl_uint N = _s -> bodies.positions.size();
    const REAL G = _s -> get_G();
    const REAL C = _s -> get_c();
    const size_t vec_size = 3 * num_bytes;

    // creates buffer for mass
    cl_mem masses_b = clCreateBuffer(context, CL_MEM_READ_ONLY, num_bytes, nullptr, &operation_result);
    Check(operation_result);

    // loads buffer for mass
    operation_result = clEnqueueWriteBuffer(queue, masses_b, CL_FALSE, 0, num_bytes, _s -> bodies.masses.data(), 0, nullptr, nullptr);
    Check(operation_result);

    // creates buffer for temperature
    cl_mem temp_b = clCreateBuffer(context, CL_MEM_READ_ONLY, num_bytes, nullptr, &operation_result);
    Check(operation_result);
    
    // loads buffer for temperature
    operation_result = clEnqueueWriteBuffer(queue, temp_b, CL_FALSE, 0, num_bytes, _s -> bodies.temps.data(), 0, nullptr, nullptr);
    Check(operation_result);

    // creates buffer for position
    cl_mem pos_b = clCreateBuffer(context, CL_MEM_READ_ONLY, vec_size, nullptr, &operation_result);
    Check(operation_result);

    // loads buffer for position
    operation_result = clEnqueueWriteBuffer(queue, pos_b, CL_FALSE, 0, vec_size, _s -> bodies.positions.data(), 0, nullptr, nullptr);
    Check(operation_result);

    // creates buffer for velocities
    cl_mem vel_b = clCreateBuffer(context, CL_MEM_READ_ONLY, vec_size, nullptr, &operation_result);
    Check(operation_result);

    // loads buffer for velocities
    operation_result = clEnqueueWriteBuffer(queue, vel_b, CL_FALSE, 0, vec_size, _s -> bodies.velocities.data(), 0, nullptr, nullptr);
    Check(operation_result);

    //creates buffer for acceleration
    cl_mem acc_br = clCreateBuffer( context, CL_MEM_WRITE_ONLY, vec_size, nullptr, &operation_result);
    Check(operation_result);

    size_t LW_size = 256;
    size_t GW_size = ((N + LW_size - 1) / LW_size) * LW_size;

    int lower = _s -> solver -> get_lower_bound(), upper = _s -> solver -> get_upper_bound();

    Check(clSetKernelArg(k, 0, sizeof(unsigned int), &N));
    Check(clSetKernelArg(k, 1, sizeof(REAL), &G));
    Check(clSetKernelArg(k, 2, sizeof(REAL), &C));
    Check(clSetKernelArg(k, 3, sizeof(cl_mem), &temp_b));
    Check(clSetKernelArg(k, 4, sizeof(cl_mem), &masses_b)); 
    Check(clSetKernelArg(k, 5, sizeof(cl_mem), &pos_b));
    Check(clSetKernelArg(k, 6, sizeof(cl_mem), &vel_b));
    Check(clSetKernelArg(k, 7, sizeof(cl_mem), &acc_br));
    Check(clSetKernelArg(k, 8, sizeof(cl_int), &lower));
    Check(clSetKernelArg(k, 9, sizeof(cl_int), &upper));

    operation_result = clEnqueueNDRangeKernel(queue, k, 1, 0, &GW_size, &LW_size, 0, nullptr, nullptr );
    Check(operation_result);

    operation_result = clEnqueueReadBuffer( queue, acc_br, CL_FALSE, 0, vec_size, _s -> bodies.accs.data(), 0, nullptr, nullptr );
    Check(operation_result);

    clFinish(queue);

    clReleaseMemObject(pos_b);
    clReleaseMemObject(acc_br);
    clReleaseMemObject(masses_b);
    clReleaseMemObject(temp_b);
    clReleaseMemObject(vel_b);
}


}
}