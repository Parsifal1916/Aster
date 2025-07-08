#pragma once 

#include <functional>

#include "Aster/physics/body.h"
#include "Aster/physics/vectors.h"


namespace Aster{
template <typename T> class BodyArray;
template <typename T> class Simulation;

template <typename T>
using func_ptr = std::function<void(Simulation<T>*)>;

template <typename T>
using force_func = T(*)(REAL, REAL, T, T, T, T, Simulation<T>*);

enum force_type:  int {NEWTON = 0, PN1 = 1, PN2 = 2, PN25 = 3, CUSTOM_F = 4};
enum update_type: int {EULER = 0, SABA2 = 1, SABA3 = 2, SABA4 = 3, SABA5 = 4, SABA6 = 5, SABA7 = 6, SABA8 = 7, SABA9 = 8, SABA10 = 9, CUSTOM_U};
enum render_style: int {DETAILED = 0, SIMPLE = 1, THERMAL = 2, TRAJ = 3};

#define SABA1 EULER
#define LEAPFROG SABA2 

enum forces_update_type: int {SINGLE_CORE, PARALLEL, GPU_UF, CUSTOM_FU};

template <typename T>
force_func<T> get_force_func(force_type type);

template <typename T>
func_ptr<T> get_update_func(update_type type, bool gpu = false);

template <typename T>
func_ptr<T> get_uf(forces_update_type type, force_type f);
}