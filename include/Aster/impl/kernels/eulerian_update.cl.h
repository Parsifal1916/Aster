#include <string>

namespace Aster{
namespace GPU{

inline std::string euler_cl = 
"#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n"
"__kernel void euler(\n"
"    const uint    N,\n"
"    const double dt,\n"
"    __global double2*  pos,\n"
"    __global double2*  vel,\n"
"    __global double2*  acc)\n"
"{\n"
"    const uint i = get_global_id(0);\n"
"    if (i >= N) return;\n"
"    vel[i] += acc[i] * dt;\n"
"    pos[i] += vel[i] * dt;\n"
"    acc[i] = (double2)(0.0,0.0);\n"
"}";

inline std::string euler_cl_lite = 
"__kernel void euler(\n"
"    const uint    N,\n"
"    const float dt,\n"
"    __global float2*  pos,\n"
"    __global float2*  vel,\n"
"    __global float2*  acc)\n"
"{\n"
"    const uint i = get_global_id(0);\n"
"    if (i >= N) return;\n"
"    vel[i] += acc[i] * dt;\n"
"    pos[i] += vel[i] * dt;\n"
"    acc[i] = (float2)(0.0,0.0);\n"
"}";

}
}
