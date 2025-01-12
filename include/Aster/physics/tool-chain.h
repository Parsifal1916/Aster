#pragma once
#include "Aster/physics/body.h"
#include "Aster/simulations/3d_sim_obj.h"
#include "Aster/simulations/sim_obj.h"

namespace Aster{

constexpr double c1 = 1 / (2 * (2 - std::pow(2, 1.0/3)));
constexpr double c2 = (1 - pow(2, 1.0/3)) * c1;
constexpr double d1 = 1 / (2 - std::pow(2, 1.0/3));
constexpr double d2 = -std::pow(2, 1.0/3) *d1;  

//===---------------------------------------------------------===//
// 2d methods                                                    //
//===---------------------------------------------------------===//

void update_euler(Body* b, Simulation* _s);
void update_leapfrog(Body* b, Simulation* _s);
void update_symplectic4(Body* body, Simulation* _s);


vec2 newtonian(double m1, double m2, vec2 v1, vec2 v2, vec2 p1, vec2 p2, Simulation* _s);

/*
* evaluates the post-newtonian approx.
* @param m1: mass of the first object
* @param m2: mass of the second object
* @param v1: velocity of the first object
* @param v2: velocity of the second object 
* @param p1: position of the first object
* @param p2: position if the second object
*/
vec2 pn2(double m1, double m2, vec2 v1, vec2 v2, vec2 p1, vec2 p2, Simulation* _s);
std::vector<vec2> get_new_pos(vec2* position, vec2* velocity, vec2* acceleration, double step);
vec2 rk4(double m1, double m2, vec2 v1, vec2 v2, vec2 p1, vec2 p2, Simulation* _s);

//===---------------------------------------------------------===//
// 3d methods                                                    //
//===---------------------------------------------------------===//

void update_euler_3d(Body3d* b, Simulation3d* _s);
void update_leapfrog_3d(Body3d* b, Simulation3d* _s);
void update_symplectic4_3d(Body3d* body,Simulation3d* _s);


vec3 newtonian_3d(double m1, double m2, vec3 v1, vec3 v2, vec3 p1, vec3 p2, Simulation3d* _s);
/*
* evaluates the post-newtonian approx.
* @param m1: mass of the first object
* @param m2: mass of the second object
* @param v1: velocity of the first object
* @param v2: velocity of the second object 
* @param p1: position of the first object
* @param p2: position if the second object
*/
vec3 pn2_3d(double m1, double m2, vec3 v1, vec3 v2, vec3 p1, vec3 p2, Simulation3d* _s);
std::vector<vec3> get_new_pos_3d(vec3* position, vec3* velocity, vec3* acceleration, double step);
vec3 rk4_3d(double m1, double m2, vec3 v1, vec3 v2, vec3 p1, vec3 p2, Simulation3d* _s);


//===---------------------------------------------------------===//
// Cosmology stuff                                               //
//===---------------------------------------------------------===//

/*
* gets the delta a from the previous 
* value of the scale factor 
*/
float get_da(float s, Simulation* _s);

/*
* updates the scale factor with the
* rk4 method
*/
void update_scale(Simulation* _s);

/*
*/


}
