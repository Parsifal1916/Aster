#include <string>

namespace Aster{
namespace GPU{

inline std::string general_saba3d = 
"#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n"
"__kernel void saba(\n"
"    const uint    N,\n"
"    const double dt,\n"
"    const double d,\n"
"    const double c,\n"
"    __global double*  pos,\n"
"    __global double*  vel,\n"
"    __global double*  acc)\n"
"{\n"
"    const uint i = get_global_id(0);\n"
"    if (i >= N) return;\n"
"    if (c != 0.0) vstore3(vload3(0, acc + 3*i) * dt * c + vload3(0, vel + i *3), 0, vel + i *3);\n"
"    vstore3(vload3(0, vel + i*3) * dt * d + vload3(0, pos + i*3), 0, pos + i *3);\n"
"    vstore3((double3)(0.0, 0.0, 0.0),0, acc + 3 *i);\n"
"}";

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