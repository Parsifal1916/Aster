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
"    __global double3*  pos,\n"
"    __global double3*  vel,\n"
"    __global double3*  acc)\n"
"{\n"
"    const uint i = get_global_id(0);\n"
"    if (i >= N) return;\n"
"    if (c != 0.0) vel[i] += acc[i] * dt * c;\n"
"    pos[i] += vel[i] * dt * d;\n"
"    acc[i] = (double3)(0.0, 0.0, 0.0);\n"
"}";


}
}