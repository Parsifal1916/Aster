#include <string>

namespace Aster{
namespace GPU{

inline std::string newton_cl = 
"#define WG 256\n"
"#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n"
"__kernel void newton(\n"
"    const uint    N,\n"
"    const double  G,\n"
"    const double  C,\n"
"    __global const double*   t,\n"
"    __global const double*   m,\n"
"    __global const double2*  pos,\n"
"    __global const double2*  vel,\n"
"    __global       double2*  acc_out)\n"
"{\n"
"    const uint i = get_global_id(0);\n"
"    if (i >= N) return;\n"
"\n"
"    double2 pi = pos[i];\n"
"    double2 ai = (double2)(0.0, 0.0);\n"
"\n"
"    __local double2 spos[WG];\n"
"    __local double  sm  [WG];\n"
"\n"
"    for (uint tile = 0; tile < N; tile += WG) {\n"
"        uint lid = get_local_id(0);\n"
"        uint j   = tile + lid;\n"
"\n"
"        if (j < N) {\n"
"            spos[lid] = pos[j];\n"
"            sm  [lid] = m[j];\n"
"        } else {\n"
"            spos[lid] = (double2)(0.0,0.0);\n"
"            sm  [lid] = 0.0;\n"
"        }\n"

"        barrier(CLK_LOCAL_MEM_FENCE);\n"
"\n"
"        for (int k = 0; k < WG; ++k) {\n"
"            uint global_j = tile + k;\n"
"            if (global_j < N && global_j != i) {\n"
"                double2 r = spos[k] - pi;\n"
"                double d2 = dot(r,r) + 10e-11;\n"
"                double invDist = native_rsqrt(d2);\n"
"                double invDist3 = invDist * invDist * invDist * G;\n" 
"                ai += sm[k] * r * invDist3;\n"
"            }\n"
"        }\n"
"        barrier(CLK_LOCAL_MEM_FENCE);\n"
"    }\n"
"\n"
"    acc_out[i] = ai;\n"
"}";

}
}
