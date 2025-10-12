#pragma once
#include <string>

namespace Aster{
namespace GPU{

std::string radix_cl = R"CLC(
#define MAX_LOCAL_SIZE  128 
#define MASK 1UL << 32

inline bool cmp(ulong a, ulong b){
    return ((a >> 32) | (a << 32)) < ((b >> 32) | (b << 32));
}

inline void swap_ulong(__global ulong *a, __global ulong *b) {
    ulong tmp = *b; *b = *a; *a = tmp;
}
inline void swap_ulocal(__local ulong *a, __local ulong *b) {
    ulong tmp = *b; *b = *a; *a = tmp;
}

inline char cmp_desc_low32(ulong a, ulong b) {
    return cmp(a, b);
}

__kernel void Sort_BitonicMergesortStart(__global const ulong* inArray,
                                         __global ulong* outArray,
                                         const uint size)
{
    __local ulong local_buffer[2 * MAX_LOCAL_SIZE];

    const uint gid = get_group_id(0);
    const uint lid = get_local_id(0);
    const uint localSize = get_local_size(0);
    
    if (localSize > MAX_LOCAL_SIZE) return;

    uint base = gid * (2u * MAX_LOCAL_SIZE);
    uint i0 = base + lid;
    uint i1 = base + lid + MAX_LOCAL_SIZE;
    
    const ulong PAD = (ulong)0; 

    local_buffer[lid] = (i0 < size) ? inArray[i0] : PAD;
    local_buffer[lid + MAX_LOCAL_SIZE] = (i1 < size) ? inArray[i1] : PAD;

    uint total = 2u * MAX_LOCAL_SIZE;
    for (uint k = 2u; k <= total; k <<= 1) {
        for (uint j = k >> 1; j > 0; j >>= 1) {
            barrier(CLK_LOCAL_MEM_FENCE);
            uint idx = 2u * lid - (lid & (j - 1u));
            if (idx + j < total) {
                char dir = ((idx & k) == 0) ? 1 : 0;
                if ( cmp(local_buffer[idx], local_buffer[idx + j])  == dir) {
                    ulong tmp = local_buffer[idx];
                    local_buffer[idx] = local_buffer[idx + j];
                    local_buffer[idx + j] = tmp;
                }
            }
        }
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    if (i0 < size) outArray[i0] = local_buffer[lid];
    if (i1 < size) outArray[i1] = local_buffer[lid + MAX_LOCAL_SIZE];
}

__kernel void Sort_BitonicMergesortGlobal(__global ulong* data,
                                          const uint size,
                                          const uint k,
                                          const uint j)
{
    uint i = get_global_id(0);
    if (i >= size) return;

    uint ixj = i ^ j;
    if (ixj <= i) return; 

    if (ixj >= size) return;

    ulong a = data[i];
    ulong b = data[ixj];

    char dir = ((i & k) == 0) ? 1 : 0; // dir==1 -> descending as implemented below

    if (cmp(a, b) == dir) {
        data[i] = b;
        data[ixj] = a;
    }
}
//////////////////////////////////////////////////////////////////

)CLC";
    
}
}