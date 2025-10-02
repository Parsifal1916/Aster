#pragma once
#include <string>

namespace Aster{
namespace GPU{

inline std::string general_saba3d = R"CLC(
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
__kernel void saba(
    const uint    N,
    const double dt,
    const double c,
    const double d,
    __global double* pos,
    __global double* vel,
    __global double* acc
){
    const uint i = get_global_id(0);
    if (i >= N) return;
    double3 v = vload3(i, vel);
    double3 a = vload3(i, acc);
    double3 p = vload3(i, pos);
    
    if (d) v += a * dt * d;
    p += v * dt * c;

    vstore3(p, i, pos);
    vstore3(v, i, vel);
    vstore3((double3)(0.0,0.0,0.0), i, acc);
})CLC";

inline std::string general_saba3d_lite = 
"__kernel void saba(\n"
"    const uint    N,\n"
"    const float dt,\n"
"    const float d,\n"
"    const float c,\n"
"    __global float3*  pos,\n"
"    __global float3*  vel,\n"
"    __global float3*  acc)\n"
"{\n"
"    const uint i = get_global_id(0);\n"
"    if (i >= N) return;\n"
"    if (c != 0.0) vel[i] += acc[i] * dt * c;\n"
"    pos[i] += vel[i] * dt * d;\n"
"    acc[i] = (float3)(0.0, 0.0, 0.0);\n"
"}";


}
}