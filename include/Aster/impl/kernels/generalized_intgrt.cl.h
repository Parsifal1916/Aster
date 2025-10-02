#pragma once
#include <string>

namespace Aster{
namespace GPU{

inline std::string general_saba = R"CLC(
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
    __kernel void saba(
        const uint    N,
        const double dt,
        const double c,
        const double d,
        __global double2* pos,
        __global double2* vel,
        __global double2* acc
    ){
        const uint i = get_global_id(0);
        if (i >= N) return;
        if (d) vel[i] += acc[i] * dt * d;
        pos[i] += vel[i] * dt * c;
        acc[i] = (double2)(0.0,0.0);
    })CLC";

inline std::string general_saba_lite = 
"__kernel void saba(\n"
"    const uint    N,\n"
"    const float dt,\n"
"    const float d,\n"
"    const float c,\n"
"    __global float2*  pos,\n"
"    __global float2*  vel,\n"
"    __global float2*  acc)\n"
"{\n"
"    const uint i = get_global_id(0);\n"
"    if (i >= N) return;\n"
"    if (c != 0) vel[i] += acc[i] * dt * c;\n"
"    pos[i] += vel[i] * dt * d;\n"
"    acc[i] = (float2)(0.0,0.0);\n"
"}";

}
}