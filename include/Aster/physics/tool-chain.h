#pragma once
#include <vector>
#include <map>
#include <string>

#include "Aster/physics/body.h"
#include "Aster/simulations/basic.h"

namespace Aster{

template <typename T> class Simulation;

extern double c1;
extern double c2;
extern double d1;
extern double d2;

//===---------------------------------------------------------===//
// 2d methods                                                    //
//===---------------------------------------------------------===//

template <typename T> class Simulation;

template <typename T>
void update_euler(Body<T>* b, Simulation<T>* _s);

template <typename T>
void update_leapfrog(Body<T>* b, Simulation<T>* _s);

template <typename T>
void update_symplectic4(Body<T>* body, Simulation<T>* _s);

template <typename T>
T newtonian(double m1, double m2, T v1, T v2, T p1, T p2, Simulation<T>* _s);

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
T pn2(double m1, double m2, T v1, T v2, T p1, T p2, Simulation<T>* _s);

template <typename T>
std::vector<T> get_new_pos(T* position, T* velocity, T* acceleration, double step);

template <typename T>
T rk4(double m1, double m2, T v1, T v2, T p1, T p2, Simulation<T>* _s);

template <typename T>
std::map<std::string, func_ptr<T>> update_funcs = {
    {"Euler", update_euler<T>} , 
    {"Leapfrog", update_leapfrog<T>},
    {"Symplectic4", update_symplectic4<T>}
};

template <typename T>
std::map<std::string, force_func<T>> force_funcs = {
    {"Newton", newtonian<T>} , 
    {"Pn2", pn2<T>}
};

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
