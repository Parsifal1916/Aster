#pragma once 

#include <cmath>
#include "Aster/physics/body.h"
#include "Aster/simulations/simulation.h"
#include <iostream>

namespace Aster{


extern const double PI;
extern const int digit_precision;
extern const double digit_coefficient;

double rng_percent();
double rng_val(double a, double b);

vec3 rng_point_in_cube(vec3 start, vec3 stop = {0,0,0});
vec2 rng_point_in_square(vec2 start, vec2 stop = {0,0});
vec2 rng_point_in_circle(double max_r, double min_r = 1);
vec3 rng_point_in_cylinder(double max_r, double min_r = 1, double thickness = 1);
vec2 rng_point_in_sphere(double max_r, double min_r = 1, double thickness = 1);

vec2 rng_vec(Simulation<vec2>* _s);
vec3 rng_vec(Simulation<vec3>* _s);


void cosmic_web(Simulation<vec2>* _s, int nums, double avr_mass);

void add_body(Simulation<vec2>* _s, double mass, vec2 pos, vec2 vel);
void add_body(Simulation<vec2>* _s, Body<vec2> b);

void add_disk(Simulation<vec2>* _s, size_t nums, vec2 center, double outer, double inner, double avr_mass = 10e6, vec2 v = {0,0});

class Composite2d* add_composite(class Composite2d* _s, size_t segments, double mass, double radius, vec2 initial_position, vec2 initial_velocity);

void add_body(Simulation<vec3>* _s, double mass, vec3 pos, vec3 vel);
void add_body(Simulation<vec3>* _s, Body<vec3> b);
vec3 rotate_point(vec3 v, double phi, double theta);
Body<vec3> make_ring_helper3d(Simulation<vec3>* _s, int outer, int inner, vec3 center, double g_pull, int thickness,  double phi, double theta, double avr_mass = 10e6, vec3 v = {0,0,0});
void add_disk(Simulation<vec3>* _s, size_t nums, vec3 center, double radius, double thickness, vec3 rotation, double avr_mass = 10e6, vec3 v = {0,0,0});
vec3 rng_3d(int max_r);
Body<vec3> pick_3d_point(Simulation<vec3>* _s, int max_r, vec3 center);
void rng_sphere(Simulation<vec3>* _s, int max_r, vec3 center, int nums);
void cosmic_web(Simulation<vec3>* _s, int nums, double avr_mass);

}//fine sim
