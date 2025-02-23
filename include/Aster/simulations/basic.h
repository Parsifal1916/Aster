#pragma once 

#include "Aster/physics/body.h"
#include "Aster/physics/vectors.h"


namespace Aster{
template <typename T> class Body;
template <typename T> class Simulation;

template <typename T>
using func_ptr = void(*)(Simulation<T>*);

template <typename T>
using force_func = T(*)(double, double, T, T, T, T, Simulation<T>*);


enum force_type:  int {NEWTON = 0, PN1 = 1, PN2 = 2, PN25 = 3};
enum update_type: int {EULER = 0, LEAPFROG = 1, SYMPLECTIC4 = 2, ADE = 3, SABA1 = 4, SABA2 = 5, SABA3 = 6, SABA4, SABA5, SABA6, SABA7, SABA8, SABA9, SABA10};

template <typename T>
force_func<T> get_force_func(force_type type);

template <typename T>
func_ptr<T> get_update_func(update_type type);
}
