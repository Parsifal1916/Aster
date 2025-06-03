#include <cmath>
#include <iostream>
#define FNL_IMPL

#include "Aster/thirdparty/FastNoiseLite.h"

#include "Aster/impl/barnes_hut.tpp"
#include "Aster/impl/BHT_impl.tpp"

#include "Aster/physics/body.h"
#include "Aster/building-api/clusters.h"


namespace Aster{
extern const double PI;
std::uniform_real_distribution<double> angle_rnd(0.0f, 360.0f);
std::uniform_real_distribution<double> normalized_rnd(0.0f, 1.f);

const int digit_precision = 4;
const double digit_coefficient = std::pow(10, digit_precision);

/*                   //===---------------------------------------------------------===//
.                    // rng generation functions                                      //
.                    //===---------------------------------------------------------===*/

/**
* @brief adds a body to the simulation
* @param _s: simulation to add the body to 
* @param mass: mass of the body
* @param pos: position of the body
* @param vel: velocity of the body
*/
void add_body(Simulation<vec2>* _s, double mass, vec2 pos, vec2 vel, double temp = 0){
    if (warn_if(!_s, "the given simulation to make add the body to is a nullptr"))
        return;

    if (warn_if(mass <= 0, "cannot have a body with negative mass"))
        return;

    _s -> bodies.positions.push_back(pos);
    _s -> bodies.velocities.push_back(vel);
    _s -> bodies.masses.push_back(mass);
    _s -> bodies.temps.push_back(temp);
    _s -> bodies.accs.push_back({0,0});
    _s -> obj ++;
}

void add_body(Simulation<vec3>* _s, double mass, vec3 pos, vec3 vel, double temp = 0){
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
vec3 rotate_point(vec3 v, double phi, double theta){
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
double rng_percent(){
    if (critical_if(!digit_coefficient, "digit coefficient has been set to zero"))
        exit(-1);
    return double(std::rand() % (int)digit_coefficient) / digit_coefficient;
}

/**
* @brief returns a random value from a to b
* @warning the order of a and b does not matter
* @param a: first value
* @param b: second value
* @return a random number 
*/
double rng_val(double a, double b){
    if (a == b) return b;

    // orders the numbers based on who's bigger
    double highest = std::max(a, b);
    double lowest = std::min(a, b);

    // gets the total range
    double range = highest - lowest;

    // both of them are negative...
    if (highest < 0){
        // ... so it swaps them
        highest = std::min(a, b);
        lowest = std::max(a,b);
        range = highest - lowest;
    }

    // calculates the random number
    double zero_2max = std::rand() % int(range * digit_coefficient);
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
vec3 rng_point_in_cube(vec3 start, vec3 stop = {0,0,0}){
    return {
        rng_val(start.x, stop.x),
        rng_val(start.y, stop.y),
        rng_val(start.z, stop.z)
    };
}

/**
* @brief returns a random point in a square with vertecies start and stop
* @param start: first vertex
* @param stop: second vertex (default = {0,0,0})
* @returns the point
*/
vec2 rng_point_in_square(vec2 start, vec2 stop = {0,0}){
    return {
        rng_val(start.x, stop.x),
        rng_val(start.y, stop.y),
    };
}

/**
* @brief returns a random point in a 2d donut
* @param max_r: outer circle
* @param min_r: inner circle
* @returns the point
*/
vec2 rng_point_in_circle(double max_r, double min_r = 1){
    //if (warn_if(max_r <= 0 || min_r <= 0 || min_r > max_r, "invalid radius value when generating random point"))
    //    return {0,0};

    double theta = rng_val(0, 2*PI);
    double r = rng_val(min_r, max_r);

    return {
        r * std::cos(theta),
        r * std::sin(theta)
    };
}

/**
* @brief returns a random point in a cylinder
* @param max_r: outer circle
* @param min_r: inner circle (default = 1)
* @param thickness: hight of the cylinder (default = 1)
*/
vec3 rng_point_in_cylinder(double max_r, double min_r = 1, double thickness = 1){
    if (warn_if(max_r <= 0 || min_r <= 0 || min_r > max_r, "invalid radius value when generating random point"))
        return {0,0,0};

    if (warn_if(thickness < 0, "negative thickness is not allowed, defaulting to its absolute value"));
        thickness = -thickness;

    // generates a random angle and radius using rng_val
    double theta = rng_val(0, 2*PI);
    double r = rng_val(min_r, max_r);

    // converts the coordinates to cartesian
    return {
        r * std::cos(theta),
        r * std::sin(theta),
        (rng_percent() - .5)* thickness // shifts the cylinder downwards
    };
}

/**
* @brief returns a random point in sphere
* @param max_r: sphere's radius 
* @param min_r: inner radius (default = 1)
*/
vec3 rng_point_in_sphere(double max_r, double min_r = 1){
    if (warn_if(max_r <= 0 || min_r <= 0 || min_r > max_r, "invalid radius value when generating random point"))
        return {0,0,0};

    double theta = rng_val(0, 2*PI);
    double phi = rng_val(0, 2*PI);
    double r = rng_val(min_r, max_r);

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
vec2 rng_vec(Simulation<vec2>* _s){
    // checks if the simulation exists
    if (warn_if(!_s, "the simulation is nullptr in \"rng_vec\"")) 
        return {0,0};

    if (critical_if(!_s -> get_width() || !_s -> get_height(), "the simulation has no space! (width = 0 or height = 0)"))
        exit(-1);

    return {
        double(std::rand() % (int)(_s  -> get_width())), 
        double(std::rand() % (int)_s  -> get_height())
    };
}

/**
* @brief picks a random point in a simulation
* @param _s: the simulation from which to pick the point
* @returns the random point
*/
vec3 rng_vec(Simulation<vec3>* _s){
    // checks if the simulation exists
    if (warn_if(!_s, "the simulation is nullptr in \"rng_vec\"")) 
        return {0,0,0};

    if (critical_if(!_s -> get_width() || !_s -> get_height() || !_s -> get_depth(), "the simulation has no space! (width = 0 or height = 0 or depth = 0)"))
        exit(-1);

    return {
        double(std::rand() % (int)_s  -> get_width()), 
        double(std::rand() % (int)_s  -> get_height()), 
        double(std::rand() % (int)_s -> get_depth())
    };
}


/*                   //===---------------------------------------------------------===//
.                    // premade clusters 2d                                           //
.                    //===---------------------------------------------------------===*/

/**
* @brief generates a disk inside a simulation
* @param _s: simulation to generate in
* @param nums: bodies to generate
* @param center: center of the cluster
* @param outer: outer radius of the disk
* @param inner: inner radius of the disk
* @param avr_mass: average mass of the bodies
* @param v: velocity of the cluster 
*/
void add_disk(Simulation<vec2>* _s, size_t nums, vec2 center, double outer, double inner, double avr_mass = 10e6, vec2 v = {0,0}){
    //if (warn_if(outer <= 0 || inner <= 0 || outer < inner, "invalid radius value when generating random point"))
    //    return;

    if (warn_if(!_s, "the simulation given to make a disk is nullptr!"))
        return;
    
    if (warn_if(avr_mass <= 0, "average mass cannot be zero"))
        return;
    
    if (warn_if(nums <= 0, "cannot generate fewer than 0 bodies!"))
        return;


    // initializes the cluster
    Cluster<vec2> cluster;
    cluster.number = nums;
    cluster.name = "Disk";

    double g_pull = nums *  _s -> get_G()  * avr_mass;

    // creates the builder lambda
    cluster.builder = [g_pull, outer, inner, v, avr_mass, center, _s ](Cluster<vec2> cl2d, size_t _) {
        vec2 pos = rng_point_in_circle(outer, inner)* _s -> get_scale(); // gets a random point inside the disk

        // generates the radius from the position
        double radius = pos.sqr_magn();

        // velocty on that point
        double magn_vel =std::sqrt(_s -> get_G() * avr_mass / radius) *100;
        vec2 vel = vec2(-pos.x/radius * magn_vel, pos.y/radius * magn_vel) + v;

        // assembles the body
        add_body(_s,
            avr_mass, 
            pos + center,
            vel
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
void cosmic_web(Simulation<vec2>* _s, int nums, double avr_mass){
    if (warn_if(!_s, "the given simulation to make a cosmic web is a nullptr"))
        return;

    if (warn_if(!nums, "cannot generate a cosmic web with fewer than zero bodies"))
        return;

    if (warn_if(avr_mass <= 0, "cosmic web cannot have less than zero mass"))
        return;

    // sets up the noise generator
    fnl_state noise = fnlCreateState();
    noise.noise_type = FNL_NOISE_OPENSIMPLEX2;

    int tot_nums = nums;

    // sets up the cluster
    Cluster<vec2> cluster;
    cluster.number = nums;
    cluster.avr_mass = avr_mass;
    cluster.size = _s -> get_center() *2;
    cluster.name = "cosmic web";


    // sets up the builder lambda
    cluster.builder = [&noise, _s, avr_mass](Cluster<vec2> cl2d, size_t _) {
        vec2 pos = rng_vec(_s); // gets a random point in the simulation

        // gets a random point until it is likely enough to generate a body
        while ((fnlGetNoise2D(&noise, pos.x, pos.y) + 1)/2 < rng_percent()) 
            pos = rng_vec(_s);

        // assembles the body
        add_body(_s, 
            avr_mass,
            pos,
            vec2({0,0}),
            (fnlGetNoise2D(&noise, pos.x, pos.y) + 1) * 5
        );
    };
    // adds the cluster to the queue
    _s -> loading_queue.add_cluster(cluster);
}

/*                   //===---------------------------------------------------------===//
.                    // premade clusters 2d                                           //
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
void add_disk(Simulation<vec3>* _s, size_t nums, vec3 center, double radius, double thickness, vec3 rotation, double avr_mass = 10e6, vec3 v = {0,0,0}){
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
    Cluster<vec3> cluster;
    cluster.number = nums;
    cluster.name = "Disk";

    double g_pull = nums *  _s -> get_G()  * avr_mass;

    // sets up the builder lambda
    cluster.builder = [g_pull, radius, v, avr_mass, center, thickness, rotation, _s ](Cluster<vec3> cl3d, size_t _) {
        vec3 pos = rng_point_in_cylinder(radius, 1, thickness); // gets a random point in the disk

        // generates the radius
        double r = pos.sqr_magn() + .1;

        // calculates the tangential velocity
        double magn_vel = g_pull*r/1000000 + std::exp(-(r*r)/80);
        vec3 vel = vec3(-pos.y/r *magn_vel, pos.x/r* magn_vel, rng_percent() - 1) + v;

        // rotates velocity and position
        pos = rotate_point(pos, rotation.x, rotation.z);
        vel = rotate_point(vel, rotation.x, rotation.z);

        // composes the body
        add_body(_s, 
            avr_mass, 
            pos + center,
            vel
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
void cosmic_web(Simulation<vec3>* _s, int nums, double avr_mass){
    // sets up the noise generator
    fnl_state noise = fnlCreateState();
    noise.noise_type = FNL_NOISE_OPENSIMPLEX2;

    int tot_nums = nums;
    
    // sets up the cluster
    Cluster<vec3> cluster;
    cluster.number = nums;
    cluster.avr_mass = avr_mass;
    cluster.size = _s -> get_center() *2; // it takes up the whole thing
    cluster.name = "cosmic web";

    // makes the builder lambda
    cluster.builder = [&noise, _s, avr_mass](Cluster<vec3> cl3d, size_t _) {
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

/**
* @brief "bakes" a simulation based on the type
* @param s: type of the simulation
* @return a pointer to the simulation
*/
Simulation<vec2>* bake(simulation_types s){
    switch (s){
     
    case LIGHT:
        return new SingleThread<vec2>();
     
    case HEAVY:
        return new Parallelized<vec2>(); 

    case BARNES_HUT:
        return new Barnes::Barnes_Hut<vec2>();

    case BH_termal:
        return new Barnes::BHT<vec2>();

    default:
        if (critical_if(true, "invalid simulation type"))
            exit(-1);    
    }
    
    return nullptr;
}

Simulation<vec3>* bake3d(simulation_types s){
    switch (s){
    
    case LIGHT:
        return new SingleThread<vec3>();
    
    case HEAVY:
        return new Parallelized<vec3>(); 

    case BARNES_HUT:
        return new Barnes::Barnes_Hut<vec3>();   

    case BH_termal:
        return new Barnes::BHT<vec3>();   

    default:
        if (critical_if(true, "invalid simulation type"))
            exit(-1);
        break;  
    }

    return nullptr;
}

}