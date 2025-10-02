#pragma once
#include <string>

namespace Aster{
namespace GPU{

inline std::string cl3d_pn25 = R"CLC(
#define WG 256
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
double get_A1(double eta, double3 v, double r_dot, double m, double r){
    return -(1 + 3*eta) * dot(v, v) + 3.0/2.0 * eta * r_dot * r_dot + 2*(2 + eta) * m / r; 
}
double get_B1(double eta){
    return 2*(2 - eta); 
}
double get_A2(double eta, double3 v, double r_dot, double m, double r){
    double retval = - eta * (3 - 4 * eta) * dot(v, v) * dot(v,v);
    retval += 1.0/2.0 * eta * (13 - 4*eta) * dot(v, v) * m / r;
    retval += 3.0/2.0 * eta * (3 - 4*eta) * dot(v, v) * r_dot * r_dot;
    retval += (2 + 25*eta + 2 *eta *eta) * r_dot * r_dot * m /r;
    retval += - 15.0/8.0 * eta * (1 - 3*eta) * r_dot*r_dot*r_dot*r_dot - 3.0/4.0 * (12 + 29 * eta) * (m /r) * (m/r);
    return retval;
}
double get_B2(double eta, double3 v, double r_dot, double m, double r){
    double retval = 0.5 * eta * (15 + 4 *eta) * dot(v, v);
    retval += -3.0/2.0 * eta * (3+ 2*eta) * r_dot *r_dot;
    retval += -1.0/2.0 * (4 + 41*eta + 8*eta*eta) * m /r;
    return retval;
}
double get_A25(double eta, double3 v, double r_dot, double m, double r){
    return 3 * dot(v,v) + 17.0/3.0 * m /r;
}
double get_B25(double eta, double3 v, double r_dot, double m, double r){
    return dot(v, v) + 3 * m /r;
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
       double3 acc = (double3)(0.0,0.0,0.0);
       double3 x = p2 - p1;
       double d2 = dot(x,x) + 10e-11;
       double invDist = native_rsqrt(d2);
       double r = sqrt(dot(x,x)) + 1e-11;
       double m = m1 + m2;
       double eta = m1 * m2 / (m*m);
       double3 v = v1 - v2;
       double3 n = x * invDist;
       double r_dot = dot(v, n);
       double invDist3 = invDist * invDist * invDist * G; 
       double a_components = get_A1(eta, v, r_dot, m, r) + get_A2(eta, v, r_dot, m, r);
       double b_components = get_B1(eta                ) + get_B2(eta, v, r_dot, m, r);
       double half_a = get_A25(eta, v, r_dot, m, r);
       double half_b = get_B25(eta, v, r_dot, m, r);
       acc += n * m2 / (r*r) * G; 
       acc += -((n * a_components + v *r_dot * b_components) * m / (r*r)) / (C*C*C) * G / m2; 
       acc += -( 8.0/5.0 * eta * (m*m) / (r*r*r) * (n * r_dot * half_a - v * half_b)) / (C*C*C*C*C*m2) * G;
       return acc;
}
__kernel void pn25(
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
    double3 acc = (double3)(0.0,0.0,0.0);
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


inline std::string cl3d_pn2 = R"CLC(
#define WG 256
#pragma OPENCL EXTENSION cl_khr_fp64 : enable

double get_A1(double eta, double3 v, double r_dot, double m, double r){
    return -(1 + 3*eta) * dot(v, v) + 3.0/2.0 * eta * r_dot * r_dot + 2*(2 + eta) * m / r; 
}

double get_B1(double eta){
    return 2*(2 - eta); 
}

double get_A2(double eta, double3 v, double r_dot, double m, double r){
    double retval = - eta * (3 - 4 * eta) * dot(v, v) * dot(v,v);
    retval += 1.0/2.0 * eta * (13 - 4*eta) * dot(v, v) * m / r;
    retval += 3.0/2.0 * eta * (3 - 4*eta) * dot(v, v) * r_dot * r_dot;
    retval += (2 + 25*eta + 2 *eta *eta) * r_dot * r_dot * m /r;
    retval += - 15.0/8.0 * eta * (1 - 3*eta) * r_dot*r_dot*r_dot*r_dot - 3.0/4.0 * (12 + 29 * eta) * (m /r) * (m/r);
    return retval;
}

double get_B2(double eta, double3 v, double r_dot, double m, double r){
    double retval = 0.5 * eta * (15 + 4 *eta) * dot(v, v);
    retval += -3.0/2.0 * eta * (3+ 2*eta) * r_dot *r_dot;
    retval += -1.0/2.0 * (4 + 41*eta + 8*eta*eta) * m /r;
    return retval;
}

double get_A25(double eta, double3 v, double r_dot, double m, double r){
    return 3 * dot(v,v) + 17.0/3.0 * m /r;
}

double get_B25(double eta, double3 v, double r_dot, double m, double r){
    return dot(v, v) + 3 * m /r;
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
       double3 acc = (double3)(0.0,0.0,0.0);
       double3 x = p2 - p1;
       double d2 = dot(x,x) + 10e-11;
       double invDist = native_rsqrt(d2);
       double r = sqrt(dot(x,x)) + 1e-11;
       double m = m1 + m2;
       double eta = m1 * m2 / (m*m);
       
       double3 v = v1 - v2;
       double3 n = x * invDist;
       double r_dot = dot(v, n);
       double invDist3 = invDist * invDist * invDist * G; 
       double a_components = get_A1(eta, v, r_dot, m, r) + get_A2(eta, v, r_dot, m, r);
       double b_components = get_B1(eta) + get_B2(eta, v, r_dot, m, r);
       acc += m2 * x * invDist3;
       acc += -((n * a_components + v * r_dot * b_components) * m / (r*r)) / (C*C*C*m1) * G;
       return acc;
}

__kernel void pn2(
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

inline std::string cl3d_pn1 = R"CLC(
#define WG 256
#pragma OPENCL EXTENSION cl_khr_fp64 : enable

double get_A1(double eta, double3 v, double r_dot, double m, double r){
    return -(1 + 3*eta) * dot(v, v) + 3/2 * eta * r_dot * r_dot + 2*(2 + eta) * m / r; 
}

double get_B1(double eta){
    return 2*(2 - eta); 
}

double get_A2(double eta, double3 v, double r_dot, double m, double r){
    double retval = - eta * (3 - 4 * eta) * dot(v, v) * dot(v,v);
    retval += 1.0/2.0 * eta * (13 - 4*eta) * dot(v, v) * m / r;
    retval += 3.0/2.0 * eta * (3 - 4*eta) * dot(v, v) * r_dot * r_dot;
    retval += (2 + 25*eta + 2 *eta *eta) * r_dot * r_dot * m /r;
    retval += - 15.0/8.0 * eta * (1 - 3*eta) * r_dot*r_dot*r_dot*r_dot - 3.0/4.0 * (12 + 29 * eta) * (m /r) * (m/r);
    return retval;
}

double get_B2(double eta, double3 v, double r_dot, double m, double r){
    double retval = 1/2 * eta * (15 + 4 *eta) * dot(v, v);
    retval += -3.0/2.0 * eta * (3+ 2*eta) * r_dot *r_dot;
    retval += -1.0/2 * (4 + 41*eta + 8*eta*eta) * m /r;
    return retval;
}

double get_A25(double eta, double3 v, double r_dot, double m, double r){
    return 3 * dot(v,v) + 17.0/3.0 * m /r;
}

double get_B25(double eta, double3 v, double r_dot, double m, double r){
    return dot(v, v) + 3 * m /r;
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
       double3 acc = (double3)(0.0,0.0,0.0);
       double3 x = p2 - p1;
       double d2 = dot(x,x) + 10e-11;
       double invDist = native_rsqrt(d2);
       double r = sqrt(dot(x,x)) + 1e-11;
       double m = m1 + m2;
       double eta = m1 * m2 / (m*m);
       
       double3 v = v1 - v2;
       double3 n = x * invDist;
       double r_dot = dot(v, n);
       double invDist3 = invDist * invDist * invDist * G; 
       double a_components = -(1 + 3*eta)*dot(v,v) + 3/2 * eta * r_dot * r_dot + 2*(2 + eta) * m * invDist;
       double b_components = get_B1(eta);
       acc += m2 * x * invDist3;
       acc += -((n * a_components + v *r_dot * b_components) * m * invDist * invDist) / (C*C*C*m1) * G;
       return acc;
}

__kernel void pn1(
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