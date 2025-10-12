#pragma once 

#include <cmath>
#include "Aster/physics/body.h"
#include <iostream>

namespace Aster{

class Simulation;
extern const REAL PI;
extern const int digit_precision;
extern const REAL digit_coefficient;

/**
* @brief rotates a point by some angles
* @param v: vector to rotate
* @param phi: z-x angle in degrees
* @param theta: y-z angle in degrees
* @return rotated vector
*/
vec3 rotate_point(vec3 v, REAL phi, REAL theta);

/**
* @brief returns a number from 0 to 1
* @return a random number from 0 to 1
*/
REAL rng_percent();

/**
* @brief returns a random value from a to b or b to a 
* @warning the order of a and b does not matter
* @param a: first value
* @param b: second value
* @return a random number 
*/
REAL rng_val(REAL a, REAL b);

/**
* @brief returns a random point in a box with vertecies start and stop
* @param start: first vertex
* @param stop: second vertex (default = {0,0,0})
* @returns the point
*/
vec3 rng_point_in_cube(vec3 start, vec3 stop = {0,0,0});

/**
* @brief returns a random point in a square with vertecies start and stop
* @param start: first vertex
* @param stop: second vertex (default = {0,0,0})
* @returns the point
*/
vec2 rng_point_in_square(vec2 start, vec2 stop = {0,0});

/**
* @brief returns a random point in a 2d donut
* @param max_r: outer circle
* @param min_r: inner circle
* @returns the point
*/
vec2 rng_point_in_circle(REAL max_r, REAL min_r = 1);

/**
* @brief returns a random point in a cylinder
* @param max_r: outer circle
* @param min_r: inner circle (default = 1)
* @param thickness: hight of the cylinder (default = 1)
*/
vec3 rng_point_in_cylinder(REAL max_r, REAL min_r = 1, REAL thickness = 1);

/**
* @brief returns a random point in sphere
* @param max_r: sphere's radius 
* @param min_r: inner radius (default = 1)
*/
vec2 rng_point_in_sphere(REAL max_r, REAL min_r = 1, REAL thickness = 1);

/**
* @brief picks a random point in a simulation
* @param _s: the simulation from which to pick the point
* @returns the random point
*/
vec3 rng_vec(Simulation* _s);

/**
* @brief generates a disk onto the given simulation
* @param _s: simulation to generate in
* @param nums: bodies to generate
* @param center: center of the disk
* @param radius: radius of the disk
* @param thickness: thickness of the disk
* @param rotation: rotation of the disk
* @param avr_mass: average mass of the bodies
* @param v: velocity of the disk
*/
void add_disk(Simulation* _s, size_t nums, vec3 center, REAL radius, REAL thickness, vec3 rotation, REAL avr_mass = 10e6, REAL center_mass = 10e20, vec3 v = {0,0,0});

/**
* @brief covers the simulation with bodies using a perlin noise
* @param _s: the simulation to cover
* @param nums: number of bodies to spawn
* @param avr_mass: average mass of each body
*/
void cosmic_web(Simulation* _s, int nums, REAL avr_mass);

/**
* @brief adds a body to the simulation
* @param _s: simulation to add the body to 
* @param mass: mass of the body
* @param pos: position of the body
* @param vel: velocity of the body
*/
void add_body(Simulation* _s, REAL mass, vec2 pos, vec2 vel, REAL temp = 0);


/**
* @brief adds a body to the simulation
* @param _s: simulation to add the body to 
* @param mass: mass of the body
* @param pos: position of the body
* @param vel: velocity of the body
*/
void add_body(Simulation* _s, REAL mass, vec3 pos, vec3 vel, REAL temp = 0);

/**
* @brief covers the simulation with bodies using a perlin noise
* @param _s: the simulation to cover
* @param nums: number of bodies to spawn
* @param avr_mass: average mass of each body
*/
void cosmic_web(Simulation* _s, int nums, REAL avr_mass);
}
