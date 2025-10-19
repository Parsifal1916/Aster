#include <string>

namespace Aster{
namespace GPU{

inline std::string bottom_up3d = R"CLC(
#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void first_pass(
    const int N,
    __global const int* left_nodes,
    __global const int* right_nodes,
    __global double* node_pos,
    __global double* n_masses,
    __global int* counters,
    __global const int* parents,
    const int start, 
    const int stop
){
    int i = get_global_id(0);
    int stop_idx = (stop < 0) ? N : stop;
    if (i >= stop_idx) return;
    if (i < start) return;
    i -= start;

    int current = i;
    while (1) {
        int parent = parents[current];
        if (parent < 0) break; 
        
        int old = atomic_inc(&counters[parent]);

        if (old != 1) {
            break;
        }

        int l = left_nodes[parent];
        int r = right_nodes[parent];

        double lm = n_masses[l];
        double rm = n_masses[r];
        double tot = lm + rm;

        n_masses[parent] = tot;

        if (tot > 0.0) {
            double3 lp = vload3(l, node_pos);
            double3 rp = vload3(r, node_pos);
            double3 res = (lp * lm + rp * rm) / tot;
            vstore3(res, parent, node_pos);
        } else {
            vstore3((double3)(0.0,0.0,0.0), parent, node_pos);
        }

        current = parent;
    }
}
)CLC";

inline std::string BHH_tree_cl3d = R"CLC(
#define MASK 1UL << 32
inline int get_up(ulong m){
    return m >> 32;
}

inline int get_down(ulong m){
    return m % (1UL << 32);
}
static inline int delta_func(uint i, uint j, uint N, __global const uint2* morton) {
    if (j >= N || i >= N) return -1;
    ulong mi = morton[i].x;
    ulong mj = morton[j].x;
    if (mi == mj) {
        return 32 + clz(morton[i].y ^ morton[j].y);
    } else {
        return clz(morton[i].x ^ morton[j].x);
    }
}

__kernel void build_tree(
    const int   N,
    __global const uint2* morton,
    __global       int*   left_nodes,
    __global       int*   right_nodes,
    __global double* node_pos,
    __global double* n_masses,
    __global const double* positions,    
    __global const double* masses,
    __global int* parents,
    __global int* counters
){
    int i = get_global_id(0);
    if (i >= N) return;

    int num_internal = N-1;
    if (i >= num_internal) return;
    counters[i] = 0;
    counters[i + N] = 0;
    int d_left  = (i == 0) ? -1 : delta_func(i, i - 1, N, morton);
    int d_right = delta_func(i, i + 1, N, morton);
    int d = (d_right > d_left) ? 1 : -1;
    int d_min = (d == 1) ? d_left : d_right;
    uint l_max = 2;
    while (1) {
        int idx = (int)i + (int)l_max * d;
        if (idx < 0 || (uint)idx >= N) break;
        if (delta_func(i, (uint)idx, N, morton) <= d_min) break;
        l_max <<= 1;
    }
    uint l = 0;
    for (uint t = l_max >> 1; t > 0; t >>= 1) {
        int test_pos = (int)i + (int)(l + t) * d;
        if (test_pos >= 0 && (uint)test_pos < N) {
            if (delta_func(i, (uint)test_pos, N, morton) > d_min) {
                l += t;
            }
        }
    }
    int j_int = (int)i + (int)l * d;
    uint j = (uint)j_int;
    int d_node = delta_func(i, j, N, morton);
    uint s = 0;
    uint range_size = (j > i) ? (j - i) : (i - j);
    uint t = (range_size + 1) / 2;
    while (t > 0) {
        int test_pos = (int)i + (int)(s + t) * d;
        if (test_pos >= 0 && (uint)test_pos < N) {
            if (delta_func(i, (uint)test_pos, N, morton) > d_node) {
                s += t;
            }
        }
        if (t == 1) break;
        t = (t + 1) / 2;
    }
    int gamma_int = (int)i + (int)s * d + ((d < 0) ? d : 0);
    uint gamma = (uint)gamma_int;
    uint internal_node = N + i;
    uint left_range  = (i < j) ? i : j;
    uint right_range = (i < j) ? j : i;
    if (left_range == gamma) {
        left_nodes[internal_node] = gamma;
    } else {
        left_nodes[internal_node] = N + gamma;
    }
    if (right_range == gamma + 1u) {
        right_nodes[internal_node] = gamma + 1u;
    } else {
        right_nodes[internal_node] = N + gamma + 1u;
    }
    const uint idx = morton[i].y;
    double3 tmp = vload3(idx, positions);
    vstore3(tmp, i, node_pos);
    n_masses[i] = masses[idx];
    if (i == N-2){
        const uint idx2 = morton[N-1].y;
        double3 tmp2 = vload3(idx2, positions);
        vstore3(tmp2, N-1, node_pos);
        n_masses[N-1] = masses[idx2];
    }
    parents[left_nodes[internal_node]] = internal_node;
    parents[right_nodes[internal_node]] = internal_node;
    parents[N] = -1;
}
)CLC";

inline std::string BHH_force_basic3d = R"CLC(
#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void barnes_force(
    const uint    N,
    const double  tree_size,
    const double  G,
    const double  C,
    const double  theta,
    __global const int*      lefts,
    __global const int*      rights,
    __global const double*   node_masses,
    __global const double*   com,
    __global const double*   pos,
    __global const double*   body_masses,
    __global       double*   acc_out,
    const int start, 
    const int stop
){
    int gid = get_global_id(0);
    const int act_N = stop - start;
    if (gid >= act_N) return;


    int max_nodes = 2* act_N -1;

    double3 my_pos = vload3(gid, pos);
    double  my_mass = body_masses[gid];
    double3 acc = (double3)(0.0, 0.0, 0.0);

    const int STACK_SIZE = 1024;
    int  stack[STACK_SIZE];
    double sizes[STACK_SIZE];
    int sp = 0;

    stack[sp] = act_N;
    sizes[sp] = tree_size;

    while (sp >= 0) {
        int node = stack[sp];
        double size = sizes[sp];
        sp--;

        if (node < 0 || node >= max_nodes) continue;

        double3 nodepos = vload3(node, com);
        double3 r = my_pos - nodepos;
        double d_squared = dot(r, r);

        bool is_leaf = (node < act_N);

        if ((size * size / d_squared < theta) || is_leaf) {
            if (d_squared > 0.0) {
                acc += get_force(G, C, my_mass, node_masses[node], my_pos, nodepos, (double3)(0.0,0.0,0.0), (double3)(0.0,0.0,0.0));
            }
            continue;
        }

        if (rights[node] != -1 && sp+1 < STACK_SIZE) {
            sp++;
            stack[sp] = rights[node];
            sizes[sp] = size * 0.5;
        }
        if (lefts[node] != -1 && sp+1 < STACK_SIZE) {
            sp++;
            stack[sp] = lefts[node];
            sizes[sp] = size * 0.5;
        }
    }

    vstore3(acc, gid, acc_out);
}
)CLC";



}}