#pragma once
#include <string>

namespace Aster{
namespace GPU{

inline std::string leapfrog_cl_3d = R"CLC(
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
__kernel void leapfrog(
    const uint    N,
    const double dt,
    const int first,
    __global double* pos,
    __global double* vel,
    __global double* acc
){
    const uint i = get_global_id(0);
    if (i >= N) return;
    double3 v = vload3(i, vel);
    double3 a = vload3(i, acc);
    double3 p = vload3(i, pos);
    
    if (first) p += v * dt + .5 * a * dt* dt;
    v += .5 * a * dt;
    
    vstore3(p, i, pos);
    vstore3(v, i, vel);
    vstore3((double3)(0.0,0.0,0.0), i, acc);
})CLC";



}
}