#pragma once
#include <string>

namespace Aster{
namespace GPU{

inline std::string cl3d_pn = R"CLC(
#define WG 256
#pragma OPENCL EXTENSION cl_khr_fp64 : enable

inline double3 get_pn25(
    const double m1, 
    const double m2, 
    const double3 v1, 
    const double3 v2, 
    const double3 p1, 
    const double3 p2,
    const double inv_r, 
    const double3 norm,
    const double G){
        
    double pseudo_newton1 = G * m1 * inv_r;
    double pseudo_newton2 = G * m2 * inv_r;

    double sq_v_diff = dot(v1 - v2, v1 - v2);

    double3 internal_1 = (v1-v2) * (-sq_v_diff + 2 * pseudo_newton1 - 8 * pseudo_newton2);

    double3 internal_2 = norm * (norm*v1 - norm*v2) * (3* sq_v_diff - 6 * pseudo_newton1 + 52.0/3.0 * pseudo_newton2);

    return (internal_1 + internal_2) * 4.0/5.0 * G *G *m1*m2 * inv_r*inv_r*inv_r ;
}

 
inline double3 get_pn2(
    const double m1, 
    const double m2, 
    const double3 v1, 
    const double3 v2, 
    const double3 p1, 
    const double3 p2,
    const double inv_r, 
    const double3 norm,
    const double G){
    const double v2_sq = dot(v2, v2);
    const double v1_sq = dot(v1, v1);

    const double v_prod = dot(v1,v2);

    double nv1 = dot(norm , v1);
    double nv2 = dot(norm , v2);
    double nv22 = dot(nv2, nv2);
    double pseudo_newtonian1 = G * m1 * inv_r;
    double pseudo_newtonian2 = G * m2 * inv_r;
    const double nv12 = dot(nv1, nv1);

    return G * m2 * inv_r * inv_r * (norm * (-2 *v2_sq * v2_sq + 4 * v2_sq * v_prod - 2 *v_prod * v_prod
    + 3.0/2.0 * v1_sq * nv22 + 9.0/2.0 * v2_sq * nv22 - 6 * v_prod * nv22
    -15.0/8.0 * nv22 * nv22 + pseudo_newtonian1 * (-15.0/4.0 * v1_sq+5.0/4.0 *v2_sq-5.0/2.0 * v_prod
    +39.0/2.0 * nv12-39.0 * nv1 * nv2+ 17.0/2.0 * nv22)
    + pseudo_newtonian2 * (4 * v2_sq- 8 * v_prod +2* nv12
    -4 * nv1 * nv2 - 6 * nv22))
    +(v1-v2) * (v1_sq * nv2+ 4 * v2_sq * nv1- 5 * v2_sq * nv2
    -4 * v_prod * nv1+4 * v_prod * nv2 - 6 * nv1 * nv22
    +9.0/2.0 * nv22 * nv2 + pseudo_newtonian1 * (-63.0/4.0 * nv1+55.0/4.0 * nv2)
    + pseudo_newtonian2 * (-2 * nv1- 2 * nv2)))
    +G*G*G * m2 * inv_r*inv_r*inv_r*inv_r * norm * (-57.0/4.0 *  m1*m1 - 9 * m2 * m2 - 69.0/2.0 * m1 * m2);
}

 
inline double3 get_pn1(
    const double m1, 
    const double m2, 
    const double3 v1, 
    const double3 v2, 
    const double3 p1, 
    const double3 p2,
    const double inv_r, 
    const double3 norm,
    const double G){

    const double v2_sq = dot(v2, v2);
    double3 a = norm * (-dot(v1, v1) - 2 * v2_sq + 4 *dot(v1,v2) + 3.0/2.0 * (dot(norm, v2) * dot(norm, v2)) + 5 * G * m1 * inv_r + 4 * G * m2 *inv_r);
    a += (v1-v2) * (4.0 * norm * v1 - 3 * dot(norm, v2));
    a = a * G * m2 * inv_r * inv_r;

    return a;
}

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
    double3 d = p2 - p1;
    double r2 = d.x*d.x + d.y*d.y + d.z*d.z + SOFTENING;
    double inv_r = 1.0/sqrt(r2);
    double inv_r2 = inv_r * inv_r;
    double3 norm = d * inv_r;

    double3 newton = norm * (G * m2 * inv_r2);
    double inv_c = 1.0 / C;

    double3 a1 = (double3)(0.0,0.0,0.0);
    double3 a2 = (double3)(0.0,0.0,0.0);
    double3 a25= (double3)(0.0,0.0,0.0);

#ifdef PN1
    a1  = inv_c * inv_c * get_pn1(m1, m2, v1, v2, p1, p2,inv_r, norm, G);
#endif
#ifdef PN2
    a2  = inv_c * inv_c * inv_c * inv_c * get_pn2(m1, m2, v1, v2, p1, p2, inv_r, norm, G);
#endif
#ifdef PN25
    a25 = inv_c * inv_c * inv_c * inv_c * inv_c * get_pn25(m1, m2, v1, v2, p1, p2, inv_r, norm, G);
#endif
    return  newton + a1 + a2 + a25;
}

__kernel void apply_forces(
        const uint N,
        const double G,
        const double C,
    __global const double*  temps,
    __global const double*  masses,
    __global const double*  positions,
    __global const double*  velocities,
    __global       double*  acc_out
){
    const uint i = get_global_id(0);
    if (i >= N) return;

    double3 acc = (double3)(0.0, 0.0, 0.0);

    double3 pi = vload3(i, positions);
    double this_mass = masses[i];
    double3 this_vel = vload3(i, velocities);

    __local double3 spos[WG];
    __local double  sm  [WG];
    __local double3 svel[WG];

    for (uint tile = 0; tile < N; tile += WG) {
        uint lid = get_local_id(0);
        uint j   = tile + lid;

        if (j < N) {
            spos[lid] = vload3(j, positions);
            sm  [lid] = masses[j];
            svel[lid] = vload3(j, velocities);
        } else {
            spos[lid] = (double3)(0.0,0.0,0.0);
            svel[lid] = (double3)(0.0,0.0,0.0);
            sm  [lid] = 0.0;
        }
        barrier(CLK_LOCAL_MEM_FENCE);

        for (int k = 0; k < WG; ++k) {
            uint global_j = tile + k;
            if (global_j >= N || global_j == i) continue;
            acc += get_force(G, C, this_mass, sm[k], pi, spos[k], this_vel, svel[k]);
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    vstore3(acc, i, acc_out);
})CLC";

}
}