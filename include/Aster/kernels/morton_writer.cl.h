#pragma once
#include <string>

namespace Aster{
namespace GPU{

std::string morton_writer = R"CLC(
#pragma OPENCL EXTENSION cl_khr_fp64 : enable

#define PRECISION_BITS 10

uint get_morton(double2 bd, double2 point) {
    point.x = (point.x + bd.x) / (2 * bd.x);
    point.y = (point.y + bd.y) / (2 * bd.y);

    point.x = max(point.x, (double)0.0);
    point.y = max(point.y, (double)0.0);

    point.x = min(point.x, (double)1.0);
    point.y = min(point.y, (double)1.0);

    uint ix = min((uint)(point.x * (1U << PRECISION_BITS)), (uint)(1u << PRECISION_BITS));
    uint iy = min((uint)(point.y * (1U << PRECISION_BITS)), (uint)(1u << PRECISION_BITS));
    
    ix = (uint)((ushort)ix);
    iy = (uint)((ushort)iy);

    ix = (ix | (ix << 8)) & 0x00FF00FF;
    ix = (ix | (ix << 4)) & 0x0F0F0F0F;
    ix = (ix | (ix << 2)) & 0x33333333;
    ix = (ix | (ix << 1)) & 0x55555555;

    iy = (iy | (iy << 8)) & 0x00FF00FF;
    iy = (iy | (iy << 4)) & 0x0F0F0F0F;
    iy = (iy | (iy << 2)) & 0x33333333;
    iy = (iy | (iy << 1)) & 0x55555555;

    return (iy << 1) | ix; 
}

__kernel void gen_mortons(
    const uint N,
    __global double2* positions,  
    __global ulong* mortons,
    __global const double* bd_arr  )   
{
    size_t gid = get_global_id(0);
    if (gid >= N) return;
    double2 p = positions[gid];

    double2 bd;
    bd.x = bd_arr[0];
    bd.y = bd_arr[1];

    mortons[gid] =  ((ulong)gid << 32)  | (get_morton(bd, p));
}
)CLC";

std::string morton_writer3d = R"CLC(
#pragma OPENCL EXTENSION cl_khr_fp64 : enable

#define PRECISION_BITS 10 

static uint expand_bits3(uint n) {
    n = (n | (n << 16)) & 0x030000FF;
    n = (n | (n <<  8)) & 0x0300F00F;
    n = (n | (n <<  4)) & 0x030C30C3;
    n = (n | (n <<  2)) & 0x09249249;
    return n;
}

ulong get_morton3(double3 bd, double3 point) {
    point.x = (point.x + bd.x) / (2.0 * bd.x);
    point.y = (point.y + bd.y) / (2.0 * bd.y);
    point.z = (point.z + bd.z) / (2.0 * bd.z);

    point.x = fmin(fmax(point.x, 0.0), 1.0);
    point.y = fmin(fmax(point.y, 0.0), 1.0);
    point.z = fmin(fmax(point.z, 0.0), 1.0);

    uint ix = (uint)(point.x * ((uint)1 << PRECISION_BITS));
    uint iy = (uint)(point.y * ((uint)1 << PRECISION_BITS));
    uint iz = (uint)(point.z * ((uint)1 << PRECISION_BITS));

    uint maxv = (uint)1 << PRECISION_BITS;
    ix = (ix < maxv ? ix : maxv);
    iy = (iy < maxv ? iy : maxv);
    iz = (iz < maxv ? iz : maxv);

    ix = (uint)((ushort) ix);
    iy = (uint)((ushort) iy);
    iz = (uint)((ushort) iz);

    uint ex = expand_bits3(ix);
    uint ey = expand_bits3(iy);
    uint ez = expand_bits3(iz);

    uint mort = (ulong)ex | ((ulong)ey << 1) | ((ulong)ez << 2);
    return mort;
}

__kernel void gen_mortons(
    const uint N,
    __global const double *positions,  
    __global ulong *mortons,
    __global const double* bd_arr,
    const int upper,
    const int lower             
) {
    size_t gid = get_global_id(0);
    if (gid >= upper || gid < lower) return; 

    double3 p = vload3(gid, positions);

    double3 bd;
    bd.x = bd_arr[0];
    bd.y = bd_arr[1];
    bd.z = bd_arr[2];

    ulong code = get_morton3(bd, p);
    mortons[gid] = ((ulong)gid << 32) | (code);
}
)CLC";

    
}
}