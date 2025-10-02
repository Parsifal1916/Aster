#pragma once
#include <string>

namespace Aster{
namespace GPU{

inline std::string newton_cl3d = R"CLC(
#define WG 256
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
double3 get_force(
   const double G,
   const double C,
   const double  m1,
   const double  m2,
   const double3 p1,
   const double3 p2,
   const double3 v1,
   const double3 v2
){
   double3 a = (double3)(0.0,0.0,0.0);
   double3 r = p2 - p1;
   double d2 = dot(r,r) + 10e-11;
   double invDist = native_rsqrt(d2);
   double invDist3 = invDist * invDist * invDist * G; 
   a += m2 * r * invDist3;
   return a;
}
__kernel void cl3d_newton(
    const uint    N,
    const double  G,
    const double  C,
    __global const double*   t,
    __global const double*   m,
    __global const double*   pos,
    __global const double*   vel,
    __global       double*   acc_out)
{
    const uint i = get_global_id(0);
    if (i >= N) return;

    double3 pi = vload3(i, pos);
    double3 ai = (double3)(0.0, 0.0, 0.0);

    __local double3 spos[WG];
    __local double  sm  [WG];

    for (uint tile = 0; tile < N; tile += WG) {
        uint lid = get_local_id(0);
        uint j   = tile + lid;

        if (j < N) {
            spos[lid] = vload3(j, pos);
            sm  [lid] = m[j];
        } else {
            spos[lid] = (double3)(0.0,0.0,0.0);
            sm  [lid] = 0.0;
        }

        barrier(CLK_LOCAL_MEM_FENCE);

        for (int k = 0; k < WG; ++k) {
            uint global_j = tile + k;
            if (global_j < N && global_j != i) {
               ai += get_force(G, C, 1.0, sm[k], pi, spos[k], vload3(0, vel), vload3(0, vel));
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    vstore3(ai, i, acc_out);
}
)CLC";

}
}
