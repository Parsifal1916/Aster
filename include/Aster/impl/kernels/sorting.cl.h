#include <string>

namespace Aster{
namespace GPU{

inline std::string sorting_kernels_cl = 
"#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable"
"#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable"
"__kernel void histogram_kernel("
"    __global const ulong* input,"
"    __global uint*        hist,"
"    const uint            shift,"
"    const uint            n"
") {"
"    uint gid = get_global_id(0);"
"    if (gid >= n) return;"
"    uint byte_val = (uint)((input[gid] >> shift) & 0xFFUL);"
"    uint group_id = get_group_id(0);"
"    atomic_inc(&hist[group_id * 256 + byte_val]);"
"}"
"__kernel void scatter_kernel("
"    __global const ulong* input,"
"    __global       ulong* output,"
"    __global const uint*  offsets_gb,"
"    const uint            shift,"
"    const uint            n"
"    uint gid = get_global_id(0);"
"    if (gid >= n) return;"
"    uint lid       = get_local_id(0);"
"    uint group_id  = get_group_id(0);"
"    uint local_size = get_local_size(0);"
"    ulong val    = input[gid];"
"    uint  byte_v = (uint)((val >> shift) & 0xFFUL);"
"    __local uint local_digits[256];"
"    local_digits[lid] = byte_v;"
"    barrier(CLK_LOCAL_MEM_FENCE);"
"    uint rank = 0;"
"    for (uint i = 0; i < lid; ++i) {"
"        rank += (local_digits[i] == byte_v) ? 1u : 0u;"
"    }"
"    barrier(CLK_LOCAL_MEM_FENCE);"
"    uint group_bucket_offset = offsets_gb[group_id * 256 + byte_v];"
"    uint pos = group_bucket_offset + rank;"
"    output[pos] = val;"
"}";
}
}