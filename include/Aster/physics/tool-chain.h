#pragma once
#include <vector>
#include <map>
#include <string>
#define CL_TARGET_OPENCL_VERSION 300
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS


#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif

#include "Aster/physics/body.h"
#include "Aster/physics/coeffs.h"
#include "Aster/simulations/basic.h"

namespace Aster{

template <typename T> class Simulation;

extern REAL c1;
extern REAL c2;
extern REAL d1;
extern REAL d2;

//===---------------------------------------------------------===//
// 2d methods                                                    //
//===---------------------------------------------------------===//

template <typename T> class Simulation;

template <typename T>
void update_euler(Simulation<T>* _s);

template <typename T>
void update_symplectic4(Simulation<T>* _s);

template <typename T>
void update_SABA1(Simulation<T>* _s);

template <typename T>
void update_SABA2(Simulation<T>* _s);

template <typename T>
void update_SABA3(Simulation<T>* _s);

template <typename T>
void update_SABA4(Simulation<T>* _s);

template <typename T>
void update_SABA4(Simulation<T>* _s);

template <typename T>
void update_SABA5(Simulation<T>* _s);

template <typename T>
void update_SABA6(Simulation<T>* _s);

template <typename T>
void update_SABA7(Simulation<T>* _s);

template <typename T>
void update_SABA8(Simulation<T>* _s);

template <typename T>
void update_SABA9(Simulation<T>* _s);

template <typename T>
void update_SABA10(Simulation<T>* _s);

template <typename T>
T newtonian(REAL m1, REAL m2, T v1, T v2, T p1, T p2, Simulation<T>* _s);

/*
* evaluates the post-newtonian approx.
* @param m1: mass of the first object
* @param m2: mass of the second object
* @param v1: velocity of the first object
* @param v2: velocity of the second object 
* @param p1: position of the first object
* @param p2: position if the second object
*/
template <typename T>
T pn2(REAL m1, REAL m2, T v1, T v2, T p1, T p2, Simulation<T>* _s);

template <typename T>
std::vector<T> get_new_pos(T* position, T* velocity, T* acceleration, REAL step);

template <typename T>
T rk4(REAL m1, REAL m2, T v1, T v2, T p1, T p2, Simulation<T>* _s);

template <typename T, typename F>
void for_each_body(Simulation<T>* _s, F func);

template <typename T>
REAL get_eccentricity(Simulation<T>* _s, size_t body, REAL relv_sq, REAL w_squared,  REAL radius, REAL mass2);

template <typename T>
void get_new_temp(Simulation<T>* _s, size_t body, size_t body2);

template <typename T>
void get_new_temp(Simulation<T>* _s, size_t body, T pos, T vel, REAL temp, REAL mass);

template <typename T>
void compute_rad_pressure(Simulation<T>* _s, size_t body, T pos, REAL temp);
 

//===---------------------------------------------------------===//
// Cosmology stuff                                               //
//===---------------------------------------------------------===//

/*
* gets the delta a from the previous 
* value of the scale factor 
*/
template <typename T>
float get_da(float s, Simulation<T>* _s);

/*
* updates the scale factor with the
* rk4 method
*/
template <typename T>
void update_scale(Simulation<T>* _s);

/*
*/

}
#include "Aster/impl/tc_impl_cl.tpp"
#include "Aster/impl/SABA.tpp"