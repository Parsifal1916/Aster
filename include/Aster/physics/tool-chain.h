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

 class Simulation;

extern REAL c1;
extern REAL c2;
extern REAL d1;
extern REAL d2;

//===---------------------------------------------------------===//
// 2d methods                                                    //
//===---------------------------------------------------------===//

 class Simulation;


void update_euler(Simulation* _s);

void update_leapfrog(Simulation* _s);

void update_WH_planetary(Simulation* _s);

void update_SABA1(Simulation* _s);
void update_SABA2(Simulation* _s);
void update_SABA3(Simulation* _s);
void update_SABA4(Simulation* _s);
void update_SABA5(Simulation* _s);
void update_SABA6(Simulation* _s);
void update_SABA7(Simulation* _s);
void update_SABA8(Simulation* _s);
void update_SABA9(Simulation* _s);
void update_SABA10(Simulation* _s);


/*
* evaluates the post-newtonian approx.
* @param m1: mass of the first object
* @param m2: mass of the second object
* @param v1: velocity of the first object
* @param v2: velocity of the second object 
* @param p1: position of the first object
* @param p2: position if the second object
*/


vec3 pn2(REAL m1, REAL m2, vec3 v1, vec3 v2, vec3 p1, vec3 p2, Simulation* _s);


vec3 pn1(REAL m1, REAL m2, vec3 v1, vec3 v2, vec3 p1, vec3 p2, Simulation* _s);


vec3 pn25(REAL m1, REAL m2, vec3 v1, vec3 v2, vec3 p1, vec3 p2, Simulation* _s);


std::vector<vec3> get_new_pos(vec3* position, vec3* velocity, vec3* acceleration, REAL step);


REAL get_eccentricity(Simulation* _s, size_t body, REAL relv_sq, REAL w_squared,  REAL radius, REAL mass2);


void get_new_temp(Simulation* _s, size_t body, vec3 pos, vec3 vel, REAL temp, REAL mass);


void compute_rad_pressure(Simulation* _s, size_t body, vec3 pos, REAL temp);
 

//===---------------------------------------------------------===//
// Cosmology stuff                                               //
//===---------------------------------------------------------===//

}