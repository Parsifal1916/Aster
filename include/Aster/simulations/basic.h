#pragma once 

#include <functional>

#include "Aster/physics/body.h"
#include "Aster/physics/vectors.h"
#include "Aster/building-api/logging.h"

#ifdef _MSC_VER
    #define FORCE_INLINE __forceinline
#else
    #define FORCE_INLINE inline __attribute__((always_inline))
#endif

#if defined(_MSC_VER)
    #define FORCE_UNROLL __pragma(loop(ivdep)) __pragma(unroll)
#elif defined(__clang__) || defined(__GNUC__)
    #define DO_PRAGMA(x) _Pragma(#x)
    #define FORCE_UNROLL DO_PRAGMA(unroll)
#else
    #define FORCE_UNROLL
#endif

namespace Aster{
#define Check(a){                       \
    if (critical_if(a != CL_SUCCESS, "error in calculating the force (GPU side)")){\
        std::cout << "error " << a << " at " << __FILE__ << ":" << __LINE__ << "\n";exit(-1);\
    }\
}

class BodyArray;
class Simulation;


using func_ptr = std::function<void(Simulation*)>;


using force_func = vec3(*)(REAL, REAL, vec3, vec3, vec3, vec3, Simulation*);

enum force_type:  int {NEWTON = 0, PN1 = 1, PN2 = 2, PN25 = 3, CUSTOM_F = 4};
enum update_type: int {EULER = 0, SABA =1, LEAPFROG=2, CUSTOM_U=3};
enum solver_type: int {SINGLE_THREAD = 0, PARALLEL = 1, BARNES_HUT = 2, MIXED_BARNES = 3, SIMPLE_GPU = 4, GPU_BARNES_HUT = 5};

}