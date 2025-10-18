#include <vector>
#include <unordered_map>
#include "Aster/simulations/basic.h"

#include "Aster/building-api/GPU_endpoint.h"

#include "Aster/simulations/sim_obj.h"
#include "Aster/simulations/barnes-hut.h"
#include "Aster/simulations/barnes_hut_GPU.h"
#include "Aster/building-api/builder.h"
#include "Aster/simulations/full_GPU.h"

#include "Aster/kernels/leapfrog.cl.h"
#include "Aster/kernels/eulerian_update.cl.h"
#include "Aster/kernels/WH_planetary.cl.h"

namespace Aster{


inline std::unordered_map<update_type, std::vector<func_ptr>> integrator_mapper = {
    {EULER, {update_euler}},
    {LEAPFROG, {update_leapfrog}},
    {WH_PLANETARY, {update_WH_planetary}},
    {SABA, {update_SABA1, update_SABA2, update_SABA3, update_SABA4, update_SABA5, update_SABA6, update_SABA7, update_SABA8, update_SABA9, update_SABA10}}
};

inline void update_euler_gpu_3d(Simulation* _s){
    _s -> solver -> compute_forces();
    using namespace GPU;
    if (!has_initialized) {
        init_opencl();
        has_initialized = true;
    }

    std::string k_name = "euler";
    static auto kernel = compile_kernel(&k_name, &eulerian_update_cl_3d);

    const size_t N = _s -> bodies.positions.size();

    size_t LW_size = 64;
    size_t GW_size = ((N + LW_size - 1) / LW_size) * LW_size;
    REAL dt = _s -> get_dt();

    Check(clSetKernelArg(kernel, 0, sizeof(unsigned int), &N));
    Check(clSetKernelArg(kernel, 1, sizeof(REAL), &dt));
    Check(clSetKernelArg(kernel, 2, sizeof(cl_mem), &_s -> positions_cl));
    Check(clSetKernelArg(kernel, 3, sizeof(cl_mem), &_s -> velocities_cl));
    Check(clSetKernelArg(kernel, 4, sizeof(cl_mem), &_s -> accs_cl));

    Check(clEnqueueNDRangeKernel(queue, kernel, 1, 0, &GW_size, &LW_size, 0, nullptr, nullptr ));
}

inline void update_leapfrog_gpu_3d(Simulation* _s){
    _s -> solver -> compute_forces();
    using namespace GPU;

    std::string k_name = "leapfrog";
    static auto kernel = compile_kernel(&k_name, &leapfrog_cl_3d);

    const size_t N = _s -> bodies.positions.size();

    size_t LW_size = 64;
    size_t GW_size = ((N + LW_size - 1) / LW_size) * LW_size;
    REAL dt = _s -> get_dt();

    int first = 1, second = 0;

    Check(clSetKernelArg(kernel, 0, sizeof(unsigned int), &N));
    Check(clSetKernelArg(kernel, 1, sizeof(REAL), &dt));
    Check(clSetKernelArg(kernel, 2, sizeof(int), &first));
    Check(clSetKernelArg(kernel, 3, sizeof(cl_mem), &_s -> positions_cl));
    Check(clSetKernelArg(kernel, 4, sizeof(cl_mem), &_s -> velocities_cl));
    Check(clSetKernelArg(kernel, 5, sizeof(cl_mem), &_s -> accs_cl));

    Check(clEnqueueNDRangeKernel(queue, kernel, 1, 0, &GW_size, &LW_size, 0, nullptr, nullptr ));

    _s -> solver -> compute_forces();

    Check(clSetKernelArg(kernel, 0, sizeof(unsigned int), &N));
    Check(clSetKernelArg(kernel, 1, sizeof(REAL), &dt));
    Check(clSetKernelArg(kernel, 2, sizeof(int), &second));
    Check(clSetKernelArg(kernel, 3, sizeof(cl_mem), &_s -> positions_cl));
    Check(clSetKernelArg(kernel, 4, sizeof(cl_mem), &_s -> velocities_cl));
    Check(clSetKernelArg(kernel, 5, sizeof(cl_mem), &_s -> accs_cl));

    Check(clEnqueueNDRangeKernel(queue, kernel, 1, 0, &GW_size, &LW_size, 0, nullptr, nullptr ));
}

inline void update_WH_planetary_gpu(Simulation* _s){
    using namespace GPU;
    std::string k_name = "wh_planetary";
    static auto kernel = compile_kernel(&k_name, &wh_cl_3d);

    const size_t N = _s -> solver -> get_range();

    size_t LW_size = 64;
    size_t GW_size = ((N + LW_size - 1) / LW_size) * LW_size;

    const REAL G = _s -> get_G(), dt = _s -> get_dt();
    const int first = 0, second = 1, lower = _s -> solver -> get_lower_bound(), upper = _s -> solver -> get_upper_bound();

    Check(clSetKernelArg(kernel, 0, sizeof(cl_mem), &_s -> positions_cl));
    Check(clSetKernelArg(kernel, 1, sizeof(cl_mem), &_s -> velocities_cl));
    Check(clSetKernelArg(kernel, 2, sizeof(cl_mem), &_s -> accs_cl));
    Check(clSetKernelArg(kernel, 3, sizeof(cl_mem), &_s -> masses_cl));
    Check(clSetKernelArg(kernel, 4, sizeof(REAL), &dt));
    Check(clSetKernelArg(kernel, 5, sizeof(REAL), &G));
    Check(clSetKernelArg(kernel, 6, sizeof(unsigned int), &lower));
    Check(clSetKernelArg(kernel, 7, sizeof(unsigned int), &upper));
    Check(clSetKernelArg(kernel, 8, sizeof(unsigned int), &first));

    Check(clEnqueueNDRangeKernel(queue, kernel, 1, 0, &GW_size, &LW_size, 0, nullptr, nullptr ));

    _s -> solver -> compute_forces();

    Check(clSetKernelArg(kernel, 0, sizeof(cl_mem), &_s -> positions_cl));
    Check(clSetKernelArg(kernel, 1, sizeof(cl_mem), &_s -> velocities_cl));
    Check(clSetKernelArg(kernel, 2, sizeof(cl_mem), &_s -> accs_cl));
    Check(clSetKernelArg(kernel, 3, sizeof(cl_mem), &_s -> masses_cl));
    Check(clSetKernelArg(kernel, 4, sizeof(REAL), &dt));
    Check(clSetKernelArg(kernel, 5, sizeof(REAL), &G));
    Check(clSetKernelArg(kernel, 6, sizeof(unsigned int), &lower));
    Check(clSetKernelArg(kernel, 7, sizeof(unsigned int), &upper));
    Check(clSetKernelArg(kernel, 8, sizeof(unsigned int), &second));

    Check(clEnqueueNDRangeKernel(queue, kernel, 1, 0, &GW_size, &LW_size, 0, nullptr, nullptr ));
}

func_ptr resolve_gpu_updater(update_type t, int ord){
    using namespace GPU;
    
    switch (t){
    case EULER:        return update_euler_gpu_3d;
    case LEAPFROG:     return update_leapfrog_gpu_3d;
    case WH_PLANETARY: return update_WH_planetary_gpu;
    case SABA:         return compile_ub_saba(ord);
    default:
        critical_if(true, "The given integrator type cannot be compiled to GPU");
        return update_euler_gpu_3d;
    }
}


func_ptr bake_update_function(update_type _t, int ord){
    int idx = static_cast<int>(_t);

    if (err_if(idx >= CUSTOM_U || idx< 0, "Invalid update type")) idx = 0; 

    if (err_if(ord >= integrator_mapper[_t].size() || ord < 0, "Cannot find the right order for the given integrator")) ord = 0;

    return integrator_mapper[_t][ord];
}


Solver* bake_solver(Simulation* _s, solver_type _t){
    int idx = static_cast<int>(_t);

    if (err_if(idx > 5 || idx< 0, "Invalid solver type")) idx = 0; 


    switch (_t) {
        case SINGLE_THREAD  : return new SingleThread(_s);
        case PARALLEL       : return new Parallelized(_s);
        case BARNES_HUT     : return new Barnes::Barnes_Hut(_s);
        case GPU_BARNES_HUT : return new Barnes::BH_hyper(_s);
        case MIXED_BARNES   : return new Barnes::BHG(_s);
        case SIMPLE_GPU     : return new SimpleGPU(_s);
    }
    return nullptr;
}

// una generica per i saba

// una per il leapfrog



}