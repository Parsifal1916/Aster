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

inline std::string eulerian_update_cl_3d = R"CLC(
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
__kernel void euler(
    const uint    N,
    const double dt,
    __global double* pos,
    __global double* vel,
    __global double* acc
){
    const uint i = get_global_id(0);
    if (i >= N) return;
    double3 v = vload3(i, vel);
    double3 a = vload3(i, acc);
    double3 p = vload3(i, pos);
    
    v += a * dt;
    p += v * dt;

    vstore3(p, i, pos);
    vstore3(v, i, vel);
    vstore3((double3)(0.0,0.0,0.0), i, acc);
})CLC";

}
}
