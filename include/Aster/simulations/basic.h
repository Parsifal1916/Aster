#pragma once 

#include "Aster/physics/body.h"
#include "Aster/physics/vectors.h"


namespace Aster{
using func_ptr3d = void(*)(class Body3d*, class Simulation3d*);
using force_func3d = vec3(*)(double, double, vec3, vec3, vec3, vec3, class Simulation3d*);

using func_ptr = void(*)(class Body*, class Simulation*);
using force_func = vec2(*)(double, double, vec2, vec2, vec2, vec2, class Simulation*);
}
