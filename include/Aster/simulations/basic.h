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
}
