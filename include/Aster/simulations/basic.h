#pragma once 

#include "Aster/physics/body.h"
#include "Aster/physics/vectors.h"


namespace Aster{
template <typename T> class Body;
template <typename T> class Simulation;

template <typename T>
using func_ptr = void(*)(class Body<T>*, Simulation<T>*);

template <typename T>
using force_func = T(*)(double, double, T, T, T, T, Simulation<T>*);


enum force_type:  int {NEWTON = 0, PN1 = 1, PN2 = 2, PN25 = 3};
enum update_type: int {EULER = 0, LEAPFROG = 1, SYMPLECTIC4 = 2};

template <typename T>
force_func<T> get_force_func(force_type type);

template <typename T>
func_ptr<T> get_update_func(update_type type);
}
