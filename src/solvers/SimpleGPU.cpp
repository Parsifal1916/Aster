#include "Aster/simulations/sim_obj.h"
#include "Aster/building-api/builder.h"
#include "Aster/building-api/GPU_endpoint.h"
#include "Aster/kernels/newton_update3d.cl.h"


namespace Aster{

SimpleGPU::SimpleGPU(Simulation* s){
    this -> _s = s;
}

void SimpleGPU::load(){
    this -> load_force_kernel = GPU::compile_uf(this -> _t, this->_s->softening);
}

void SimpleGPU::compute_forces(){
    using namespace GPU;
    static int a = 0;
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
    const size_t vec_size = N * sizeof(vec3);

    std::string name = "cl3d_newton";
    static cl_kernel k = compile_kernel(&name, &newton_cl3d, this -> _s -> softening, false);

    // creates buffer for temperature
    cl_mem temp_b = clCreateBuffer(context, CL_MEM_READ_ONLY, num_bytes, nullptr, &operation_result);
    Check(operation_result);
    
    // loads buffer for temperature
    operation_result = clEnqueueWriteBuffer(queue, temp_b, CL_FALSE, 0, num_bytes, _s -> bodies.temps.data(), 0, nullptr, nullptr);
    Check(operation_result);


    size_t LW_size = 256;
    size_t GW_size = ((N + LW_size - 1) / LW_size) * LW_size;

    int lower = get_lower_bound(), upper = get_upper_bound();

    Check(clSetKernelArg(k, 0, sizeof(unsigned int), &N));
    Check(clSetKernelArg(k, 1, sizeof(REAL), &G));
    Check(clSetKernelArg(k, 2, sizeof(REAL), &C));
    Check(clSetKernelArg(k, 3, sizeof(cl_mem), &temp_b));
    Check(clSetKernelArg(k, 4, sizeof(cl_mem), &this -> _s -> masses_cl)); 
    Check(clSetKernelArg(k, 5, sizeof(cl_mem), &this -> _s -> positions_cl));
    Check(clSetKernelArg(k, 6, sizeof(cl_mem), &this -> _s -> velocities_cl));
    Check(clSetKernelArg(k, 7, sizeof(cl_mem), &this -> _s -> accs_cl));
    Check(clSetKernelArg(k, 8, sizeof(cl_int), &lower));
    Check(clSetKernelArg(k, 9, sizeof(cl_int), &upper));

    operation_result = clEnqueueNDRangeKernel(queue, k, 1, 0, &GW_size, &LW_size, 0, nullptr, nullptr );
    Check(operation_result);

    clReleaseMemObject(temp_b);
    //Check(clEnqueueReadBuffer(queue, this -> _s -> accs_cl, CL_TRUE, 0, N*sizeof(vec3), this -> _s -> bodies.accs.data(), 0, nullptr, nullptr));
//
    //for (const auto& v : this -> _s -> bodies.accs){
    //    std::cout << v << "\n";
    //}
    //if (a > -1) exit(0);
}



}