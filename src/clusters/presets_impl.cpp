#include <cmath>
#include <iostream>
#include <random>
#define FNL_IMPL

#include "Aster/building-api/presets.h"
#include "Aster/building-api/logging.h"
#include "Aster/simulations/sim_obj.h"
#include "Aster/thirdparty/FastNoiseLite.h"

namespace Aster{
extern const REAL PI;
std::uniform_real_distribution<REAL> angle_rnd(0.0f, 360.0f);
std::uniform_real_distribution<REAL> normalized_rnd(0.0f, 1.f);

const int digit_precision = 4;
const REAL digit_coefficient = std::pow(10, digit_precision);

/*                   //===---------------------------------------------------------===//
.                    // rng generation functions                                      //
.                    //===---------------------------------------------------------===*/


void add_body(Simulation* _s, REAL mass, vec3 pos, vec3 vel, REAL temp){
    if (warn_if(!_s, "the given simulation to make add the body to is a nullptr"))
        return;

    if (warn_if(mass <= 0, "cannot have a body with negative mass"))
        return;

    _s -> bodies.positions.push_back(pos);
    _s -> bodies.velocities.push_back(vel);
    _s -> bodies.masses.push_back(mass);
    _s -> bodies.temps.push_back(temp);
    _s -> bodies.accs.push_back({0,0,0});
}



/**
* @brief rotates a point by some angles
* @param v: vector to rotate
* @param phi: z-x angle in degrees
* @param theta: y-z angle in degrees
* @return rotated vector
*/
vec3 rotate_point(vec3 v, REAL phi, REAL theta){
    float x1, z1, y1, z2;

    // converts the angles in radiants
    theta *= PI/180 ;
    phi   *= PI/180 ; 

    // applyies the first rotation
    x1 = v.x * std::cos(phi) + v.z * std::sin(phi);
    z1 = - v.x * std::sin(phi) + v.z *std::cos(phi);

    // second rotation
    y1 = v.y * std::cos(theta) + z1 *std::sin(theta);
    z2 = -v.y *std::sin(theta) + z1 *std::cos(theta);
    
    // returns the new vector
    return {x1, y1, z2};
}

/**
* @brief returns a number from 0 to 1
* @return a random number from 0 to 1
*/
REAL rng_percent(){
    if (critical_if(!digit_coefficient, "digit coefficient has been set to zero"))
        exit(-1);
    return REAL(std::rand() % (int)digit_coefficient) / digit_coefficient;
}

/**
* @brief returns a random value from a to b
* @warning the order of a and b does not matter
* @param a: first value
* @param b: second value
* @return a random number 
*/
REAL rng_val(REAL a, REAL b){
    if (a == b) return b;

    // orders the numbers based on who's bigger
    REAL highest = std::max(a, b);
    REAL lowest = std::min(a, b);

    // gets the total range
    REAL range = highest - lowest;

    // both of them are negative...
    if (highest < 0){
        // ... so it swaps them
        highest = std::min(a, b);
        lowest = std::max(a,b);
        range = highest - lowest;
    }

    // calculates the random number
    REAL zero_2max = std::rand() % int(range * digit_coefficient);
    zero_2max *= (range < 0 && lowest < 0) ? -1 : 1;

    // it then gets adjusted and scaled
    zero_2max /= digit_coefficient;
    return zero_2max + lowest;
}

/**
* @brief returns a random point in a box with vertecies start and stop
* @param start: first vertex
* @param stop: second vertex (default = {0,0,0})
* @returns the point
*/
vec3 rng_point_in_cube(vec3 start, vec3 stop){
    return {
        rng_val(start.x, stop.x),
        rng_val(start.y, stop.y),
        rng_val(start.z, stop.z)
    };
}

/**
* @brief returns a random point in a cylinder
* @param max_r: outer circle
* @param min_r: inner circle (default = 1)
* @param thickness: hight of the cylinder (default = 1)
*/
vec3 rng_point_in_cylinder(REAL max_r, REAL min_r, REAL thickness){
    if (warn_if(max_r <= 0 || min_r <= 0 || min_r > max_r, "invalid radius value when generating random point"))
        return {0,0,0};

    if (warn_if(thickness < 0, "negative thickness is not allowed, defaulting to its absolute value"));
        thickness = -thickness;

    // generates a random angle and radius using rng_val
    REAL theta = rng_val(0, 2*PI);
    REAL r = rng_val(min_r, max_r);

    // converts the coordinates to cartesian
    return {
        r * std::cos(theta),
        r * std::sin(theta),
        (rng_percent() - (REAL).5)* thickness // shifts the cylinder downwards
    };
}

/**
* @brief returns a random point in sphere
* @param max_r: sphere's radius 
* @param min_r: inner radius (default = 1)
*/
vec3 rng_point_in_sphere(REAL max_r, REAL min_r){
    if (warn_if(max_r <= 0 || min_r <= 0 || min_r > max_r, "invalid radius value when generating random point"))
        return {0,0,0};

    REAL theta = rng_val(0, 2*PI);
    REAL phi = rng_val(0, 2*PI);
    REAL r = rng_val(min_r, max_r);

    return {
        r * std::sin(theta) * std::cos(phi),
        r * std::sin(theta) * std::sin(phi),
        r * std::cos(theta)
    };
}



/**
* @brief picks a random point in a simulation
* @param _s: the simulation from which to pick the point
* @returns the random point
*/
vec3 rng_vec(Simulation* _s){
    // checks if the simulation exists
    if (warn_if(!_s, "the simulation is nullptr in \"rng_vec\"")) 
        return {0,0,0};

    if (critical_if(!_s -> get_width() || !_s -> get_height() || !_s -> get_depth(), "the simulation has no space! (width = 0 or height = 0 or depth = 0)"))
        exit(-1);

    return {
        REAL(rng_percent() *  (int)_s  -> get_width()) , 
        REAL(rng_percent() *  (int)_s  -> get_height()), 
        REAL(rng_percent() *  (int)_s -> get_depth())  
    };
}


/*                   //===---------------------------------------------------------===//
.                    // premade clusters 3d                                           //
.                    //===---------------------------------------------------------===*/

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
void add_disk(Simulation* _s, size_t nums, vec3 center, REAL radius, REAL thickness, vec3 rotation, REAL avr_mass, REAL center_mass, vec3 v){
    if (warn_if(!_s, "the given simulation to make the disk is a nullptr"))
        return;

    if (warn_if(!nums, "cannot generate a disk with fewer than zero bodies"))
        return;

    if (warn_if(avr_mass <= 0, "disk cannot have less than zero mass"))
        return;

    if (warn_if(radius <= 0, "disk cannot have negative radius"))
        return;

    if (warn_if(thickness <= 0, "disk cannot have negative thickness"))
        return;

    // sets up the cluster
    Cluster cluster;
    cluster.number = nums;
    cluster.name = "Disk";
    add_body(_s, 
        center_mass, 
        rotate_point(center, rotation.x, rotation.z),
        rotate_point(v, rotation.x, rotation.z)
    );

    // sets up the builder lambda
    cluster.builder = [radius, v, avr_mass, center, thickness, rotation, _s, center_mass ](Cluster cl3d, size_t _) {
        vec3 pos = rng_point_in_cylinder(radius, 10, thickness)* _s -> get_scale(); // gets a random point inside the disk
        //pos.y /= 2.5;
        // generates the radius from the position
        REAL r = pos.magnitude();

        // velocty on that point
        REAL magn_vel = std::sqrt(_s -> get_G() * center_mass/ r);
        vec3 vel = vec3(pos.y/r * magn_vel, -pos.x/r * magn_vel, 0);
       

        // assembles the body
        add_body(_s,
            rng_percent() * avr_mass, 
            rotate_point(pos + center, rotation.x, rotation.z),
            rotate_point(vel + v, rotation.x, rotation.z)
        );
    };

    // adds the cluster to the queue
    _s -> loading_queue.add_cluster(cluster);
}

/**
* @brief covers the simulation with bodies using a perlin noise
* @param _s: the simulation to cover
* @param nums: number of bodies to spawn
* @param avr_mass: average mass of each body
*/
void cosmic_web(Simulation* _s, int nums, REAL avr_mass){
    // sets up the noise generator
    fnl_state noise = fnlCreateState();
    noise.noise_type = FNL_NOISE_OPENSIMPLEX2;

    int tot_nums = nums;
    
    // sets up the cluster
    Cluster cluster;
    cluster.number = nums;
    cluster.avr_mass = avr_mass;
    cluster.size = _s -> get_center() *2; // it takes up the whole thing
    cluster.name = "cosmic web";

    // makes the builder lambda
    cluster.builder = [&noise, _s, avr_mass](Cluster cl3d, size_t _) {
        vec3 pos = rng_vec(_s);

        // generates a random point until it can instatiate a body
        while ((fnlGetNoise3D(&noise, pos.x, pos.y, pos.z) + 1)/2 < rng_percent()) 
            pos = rng_vec(_s);

        // assembles a body
        add_body(_s, 
            avr_mass,
            pos,
            vec3({0,0,0}),
            (fnlGetNoise3D(&noise, pos.x, pos.y, pos.z) + 1) * 5
        );
    };

    // adds it to the queue
    _s -> loading_queue.add_cluster(cluster);
}

//===---------------------------------------------------------===//
// baking                                                        //
//===---------------------------------------------------------===//

}