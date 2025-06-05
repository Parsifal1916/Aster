#include <string>

namespace Aster{
namespace GPU{

inline std::string sorting_kernels_cl = 
"#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable\n"
"#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable\n"
"__kernel void histogram_kernel(\n"
"    __global const ulong* input,\n"
"    __global uint*        hist,\n"
"    const uint            shift,\n"
"    const uint            n\n"
") {\n"
"    uint gid = get_global_id(0);\n"
"    if (gid >= n) return;\n"
"    uint byte_val = (uint)((input[gid] >> shift) & 0xFFUL);\n"
"    uint group_id = get_group_id(0);\n"
"    atomic_inc(&hist[group_id * 256 + byte_val]);\n"
"}\n"
"\n"
"__kernel void scatter_kernel(\n"
"    __global const ulong* input,\n"
"    __global       ulong* output,\n"
"    __global const uint*  offsets_gb,\n"
"    const uint            shift,\n"
"    const uint            n\n"
") {\n"
"    uint gid = get_global_id(0);\n"
"    if (gid >= n) return;\n"
"    uint lid        = get_local_id(0);\n"
"    uint group_id   = get_group_id(0);\n"
"    ulong val       = input[gid];\n"
"    uint byte_val   = (uint)((val >> shift) & 0xFFUL);\n"
"\n"
"    __local uint local_digits[256];\n"
"    local_digits[lid] = byte_val;\n"
"    barrier(CLK_LOCAL_MEM_FENCE);\n"
"\n"
"    uint rank = 0;\n"
"    for (uint i = 0; i < lid; ++i) {\n"
"        rank += (local_digits[i] == byte_val) ? 1u : 0u;\n"
"    }\n"
"    barrier(CLK_LOCAL_MEM_FENCE);\n"
"\n"
"    uint group_bucket_offset = offsets_gb[group_id * 256 + byte_val];\n"
"    uint pos = group_bucket_offset + rank;\n"
"    output[pos] = val;\n"
"}\n";

}
}