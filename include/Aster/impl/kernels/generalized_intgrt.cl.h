#include <string>

namespace Aster{
namespace GPU{

inline std::string general_saba = 
"#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n"
"__kernel void saba(\n"
"    const uint    N,\n"
"    const double dt,\n"
"    const double d,\n"
"    const double c,\n"
"    __global double2*  pos,\n"
"    __global double2*  vel,\n"
"    __global double2*  acc)\n"
"{\n"
"    const uint i = get_global_id(0);\n"
"    if (i >= N) return;\n"
"    if (c != 0) vel[i] += acc[i] * dt * c;\n"
"    pos[i] += vel[i] * dt * d;\n"
"    acc[i] = (double2)(0.0,0.0);\n"
"}";

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