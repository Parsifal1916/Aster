#include <string>

namespace Aster{
namespace GPU{

inline std::string cl_pns = 
"#define WG 256\n"
"double get_A1(double eta, double2 v, double r_dot, double m, double r){\n"
"    return -(1 + 3*eta) * dot(v, v) + 3/2 * eta * r_dot * r_dot + 2*(2 + eta) * m / r; \n"
"}\n"
"\n"
"double get_B1(double eta){\n"
"    // 2(2 - η)\n"
"    return 2*(2 - eta); \n"
"}\n"
"\n"
"double get_A2(double eta, double2 v, double r_dot, double m, double r){\n"
"    double retval = - eta * (3 - 4 * eta) * dot(v, v) * dot(v,v);\n"
"    retval += 1.0/2.0 * eta * (13 - 4*eta) * dot(v, v) * m / r;\n"
"    retval += 3.0/2.0 * eta * (3 - 4*eta) * dot(v, v) * r_dot * r_dot;\n"
"    retval += (2 + 25*eta + 2 *eta *eta) * r_dot * r_dot * m /r;\n"
"    retval += - 15.0/8.0 * eta * (1 - 3*eta) * r_dot*r_dot*r_dot*r_dot - 3.0/4.0 * (12 + 29 * eta) * (m /r) * (m/r);\n"
"\n"
"    return retval;\n"
"}\n"
"\n"
"double get_B2(double eta, double2 v, double r_dot, double m, double r){\n"
"    double retval = 1/2 * eta * (15 + 4 *eta) * dot(v, v);\n"
"    retval += -3.0/2.0 * eta * (3+ 2*eta) * r_dot *r_dot;\n"
"    retval += -1.0/2 * (4 + 41*eta + 8*eta*eta) * m /r;\n"
"\n"
"    return retval;\n"
"}\n"
"\n"
"double get_A25(\n"
"    double eta, \n"
"    double2 v, \n"
"    double r_dot, \n"
"    double m, \n"
"    double r\n"
"){\n"
"    return 3 * dot(v,v) + 17.0/3.0 * m /r;\n"
"}\n"
"\n"
"double get_B25(double eta, double2 v, double r_dot, double m, double r){\n"
"    return dot(v, v) + 3 * m /r;\n"
"}\n"
"\n"
"__kernel void pn25(\n"
"             const uint N,\n"
"             const double G,\n"
"             const double C,\n"
"    __global const double*  temps,\n"
"    __global const double*  masses,\n"
"    __global const double2* positions,\n"
"    __global const double2* velocities,\n"
"    __global       double2* acc_out\n"
"){\n"
"\n"
"    const uint i = get_global_id(0);\n"
"    if (i >= N) return;\n"
"\n"
"\n"
"    double2 acc = (double2)(0.0, 0.0);\n"
"\n"
"    double2 pi = positions[i];\n"
"    double this_mass = masses[i];\n"
"    double2 this_vel = velocities[i];\n"
"\n"
"    __local double2 spos[WG];\n"
"    __local double  sm  [WG];\n"
"    __local double2  svel[WG];\n"
"\n"
"    for (uint tile = 0; tile < N; tile += WG) {\n"
"        uint lid = get_local_id(0);\n"
"        uint j   = tile + lid;\n"
"\n"
"        if (j < N) {\n"
"            spos[lid] = positions[j];\n"
"            sm  [lid] = masses[j];\n"
"            svel[lid] = velocities[i];\n"
"        } else {\n"
"            spos[lid] = (double2)(0.0,0.0);\n"
"            svel[lid] = (double2)(0.0,0.0);\n"
"            sm  [lid] = 0.0;\n"
"        }\n"
"        barrier(CLK_LOCAL_MEM_FENCE);\n"
"\n"
"        for (int k = 0; k < WG; ++k) {\n"
"            uint global_j = tile + k;\n"
"            if (global_j >= N || global_j == i) continue;\n"
"            double2 x = (spos[k] - pi + 1e-11);\n"
"            double r = sqrt(dot(x,x)) + 1e-11;\n"
"            double m = sm[k] + this_mass;\n"
"            double eta = sm[k] * this_mass / (m*m);\n"
"            \n"
"            double2 v = this_vel - svel[k];\n"
"            double2 n = x / r;\n"
"            double r_dot = dot(v, n);\n"
"\n"
"            double a_components = get_A1(eta, v, r_dot, m, r) + get_A2(eta, v, r_dot, m, r);\n"
"            double b_components = get_B1(eta                ) + get_B2(eta, v, r_dot, m, r);\n"
"\n"
"            double half_a = get_A25(eta, v, r_dot, m, r);\n"
"            double half_b = get_B25(eta, v, r_dot, m, r);\n"
"\n"
"            acc += n * sm[k] / (r*r) * G; \n"
"            acc += -((n * a_components + v *r_dot * b_components) * m / (r*r)) / (C*C*C) * G / sm[k]; \n"
"            acc += -( 8.0/5.0 * eta * (m*m) / (r*r*r) * (n * r_dot * half_a - v * half_b)) / (C*C*C*C*C*sm[k]) * G;\n"
"            \n"
"        }\n"
"        barrier(CLK_LOCAL_MEM_FENCE);\n"
"    }\n"
"\n"
"    acc_out[i] = acc;\n"
"}\n"
"\n"
"__kernel void pn2(\n"
"             const uint N,\n"
"             const double G,\n"
"             const double C,\n"
"    __global const double*  temps,\n"
"    __global const double*  masses,\n"
"    __global const double2* positions,\n"
"    __global const double2* velocities,\n"
"    __global       double2* acc_out\n"
"){\n"
"\n"
"    const uint i = get_global_id(0);\n"
"    if (i >= N) return;\n"
"\n"
"\n"
"    double2 acc = (double2)(0.0, 0.0);\n"
"\n"
"    double2 pi = positions[i];\n"
"    double this_mass = masses[i];\n"
"    double2 this_vel = velocities[i];\n"
"\n"
"    __local double2 spos[WG];\n"
"    __local double  sm  [WG];\n"
"    __local double2  svel[WG];\n"
"\n"
"    for (uint tile = 0; tile < N; tile += WG) {\n"
"        uint lid = get_local_id(0);\n"
"        uint j   = tile + lid;\n"
"\n"
"        if (j < N) {\n"
"            spos[lid] = positions[j];\n"
"            sm  [lid] = masses[j];\n"
"            svel[lid] = velocities[j];\n"
"        } else {\n"
"            spos[lid] = (double2)(0.0,0.0);\n"
"            svel[lid] = (double2)(0.0,0.0);\n"
"            sm  [lid] = 0.0;\n"
"        }\n"
"        barrier(CLK_LOCAL_MEM_FENCE);\n"
"\n"
"        for (int k = 0; k < WG; ++k) {\n"
"            uint global_j = tile + k;\n"
"            if (global_j >= N || global_j == i) continue;\n"
"            double2 x = (spos[k] - pi);\n"
"            double r = sqrt(dot(x, x)) + 1e-11;\n"
"            double m = sm[k] + this_mass;\n"
"            double eta = this_mass*sm[k] / (m*m);\n"
"            \n"
"            double2 v = this_vel - svel[k];\n"
"            double2 n = x / r;\n"
"            double r_dot = dot(v, n);\n"
"        \n"
"            double a_components = get_A1(eta, v, r_dot, m, r) + get_A2(eta, v, r_dot, m, r);\n"
"            double b_components = get_B1(eta                ) + get_B2(eta, v, r_dot, m, r);\n"
"        \n"
"        \n"
"            acc += n * sm[k] / (r*r) *G; \n"
"            acc +=  -((n * a_components + v *r_dot * b_components) * m / (r*r)) / (C*C*C*this_mass) * G;  \n"
"            \n"
"        }\n"
"        barrier(CLK_LOCAL_MEM_FENCE);\n"
"    }\n"
"\n"
"    acc_out[i] = acc;\n"
"}\n"
"\n"
"\n"
"__kernel void pn1(\n"
"        const uint N,\n"
"        const double G,\n"
"        const double C,\n"
"    __global const double*  temps,\n"
"    __global const double*  masses,\n"
"    __global const double2* positions,\n"
"    __global const double2* velocities,\n"
"    __global       double2* acc_out\n"
"){\n"
"\n"
"    const uint i = get_global_id(0);\n"
"    if (i >= N) return;\n"
"\n"
"\n"
"    double2 acc = (double2)(0.0, 0.0);\n"
"\n"
"    double2 pi = positions[i];\n"
"    double this_mass = masses[i];\n"
"    double2 this_vel = velocities[i];\n"
"\n"
"    __local double2 spos[WG];\n"
"    __local double  sm  [WG];\n"
"    __local double2 svel[WG];\n"
"\n"
"    for (uint tile = 0; tile < N; tile += WG) {\n"
"        uint lid = get_local_id(0);\n"
"        uint j   = tile + lid;\n"
"\n"
"        if (j < N) {\n"
"            spos[lid] = positions[j];\n"
"            sm  [lid] = masses[j];\n"
"            svel[lid] = velocities[j];\n"
"        } else {\n"
"            spos[lid] = (double2)(0.0,0.0);\n"
"            svel[lid] = (double2)(0.0,0.0);\n"
"            sm  [lid] = 0.0;\n"
"        }\n"
"        barrier(CLK_LOCAL_MEM_FENCE);\n"
"\n"
"        for (int k = 0; k < WG; ++k) {\n"
"            uint global_j = tile + k;\n"
"            if (global_j >= N || global_j == i) continue;\n"
"            double2 x = spos[k] - pi;\n"
"            double r = sqrt(dot(x, x)) + 1e-11;\n"
"            double m = sm[k] + this_mass;\n"
"            double eta = this_mass*sm[k] / (m*m);\n"
"            \n"
"            double2 v = this_vel - svel[k];\n"
"            double2 n = x / r;\n"
"            double r_dot = dot(v,n);"
"        \n"
"            double a_components = -(1 + 3*eta)*dot(v,v) + 3/2 * eta * r_dot * r_dot + 2*(2 + eta) * m / r;\n"
"            double b_components = get_B1(eta                );\n"
"        \n"
"        \n"
"            acc += n * sm[k] / (r*r) * G; \n"
"            acc += -((n * a_components + v *r_dot * b_components) * m / (r*r)) / (C*C*C*this_mass) * G;// \n"
"        \n"
"        }\n"
"        barrier(CLK_LOCAL_MEM_FENCE);\n"
"        }\n"
"\n"
"    acc_out[i] = acc;\n"
"};";

//def a(text): 
//    splitted = text.split("\n")
//    for i in range(len(splitted)):
//        splitted[i] = "\"" + splitted[i] + "\\n\""
//    return "\n".join(splitted)

}
}