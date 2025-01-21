#include <cmath>
#include <iostream>
#include <noise/noise.h>

#include "Aster/physics/body.h"
#include "Aster/simulations/simulation.h"
#include "Aster/building-api/clusters.h"

namespace Aster{

noise::module::Perlin perlin;

const double PI = 3.141592653589793238462643383279;
const int digit_precision = 4;
const double digit_coefficient = std::pow(10, digit_precision);

/*                   //===---------------------------------------------------------===//
.                    // rng generation functions                                      //
.                    //===---------------------------------------------------------===*/

vec3 rotate_point(vec3 v, double phi, double theta){
    float x1, z1, y1, z2;
    theta *= PI/180 ;
    phi   *= PI/180 ; 

    x1 = v.x * std::cos(phi) + v.z * std::sin(phi);
    z1 = - v.x * std::sin(phi) + v.z *std::cos(phi);

    y1 = v.y * std::cos(theta) + z1 *std::sin(theta);
    z2 = -v.y *std::sin(theta) + z1 *std::cos(theta);

    return {x1, y1, z2};
}

double rng_percent(){
    return double(std::rand() % (int)digit_coefficient) / digit_coefficient;
}

double rng_val(double a, double b){
    if (a == b) return b;

    double highest = std::max(a, b);
    double lowest = std::min(a, b);

    double range = highest - lowest;

    if (highest < 0){
        highest = std::min(a, b);
        lowest = std::max(a,b);
        range = highest - lowest;
    }

    double zero_2max = std::rand() % int(range * digit_coefficient);
    zero_2max *= (range < 0 && lowest < 0) ? -1 : 1;

    zero_2max /= digit_coefficient;
    return zero_2max + lowest;
}

vec3 rng_point_in_cube(vec3 start, vec3 stop = {0,0,0}){
    return {
        rng_val(start.x, stop.x),
        rng_val(start.y, stop.y),
        rng_val(start.z, stop.z)
    };
}

vec2 rng_point_in_square(vec2 start, vec2 stop = {0,0}){
    return {
        rng_val(start.x, stop.x),
        rng_val(start.y, stop.y),
    };
}

vec2 rng_point_in_circle(double max_r, double min_r = 1){
    double theta = rng_val(0, 2*PI);
    double r = rng_val(min_r, max_r);

    return {
        r * std::cos(theta),
        r * std::sin(theta)
    };
}

vec3 rng_point_in_cylinder(double max_r, double min_r = 1, double thickness = 1){
    double theta = rng_val(0, 2*PI);
    double r = rng_val(min_r, max_r);

    return {
        r * std::cos(theta),
        r * std::sin(theta),
        (rng_percent() - .5)* thickness
    };
}

vec3 rng_point_in_sphere(double max_r, double min_r = 1){
    double theta = rng_val(0, 2*PI);
    double phi = rng_val(0, 2*PI);
    double r = rng_val(min_r, max_r);

    return {
        r * std::sin(theta) * std::cos(phi),
        r * std::sin(theta) * std::sin(phi),
        r * std::cos(theta)
    };
}

vec2 rng_vec2(Simulation* _s){
    return {
        double(std::rand() % (int)(_s -> data.WIDTH)), 
        double(std::rand() % (int)_s -> data.HEIGHT)
    };
}

vec3 rng_vec3(Simulation3d* _s){
    return {
        double(std::rand() % (int)_s -> data.WIDTH), 
        double(std::rand() % (int)_s -> data.HEIGHT), 
        double(std::rand() % (int)_s -> data.depth)
    };
}


/*                   //===---------------------------------------------------------===//
.                    // premade clusters 2d                                           //
.                    //===---------------------------------------------------------===*/

void add_disk(Simulation* _s, size_t nums, vec2 center, double outer, double inner, double avr_mass = 10e6, vec2 v = {0,0}){
    _s -> bodies.reserve(_s -> bodies.size() + nums);

    Cluster2d cluster;
    cluster.number = nums;
    cluster.name = "Disk";

    double g_pull = nums *  _s -> data.G  * avr_mass;

    cluster.builder = [g_pull, outer, inner, v, avr_mass, center ](Cluster2d cl2d, size_t _) {
        vec2 pos = rng_point_in_circle(outer, inner);

        double radius = pos.sqr_magn();
        double magn_vel = g_pull*radius/1000 + std::exp(-(radius*radius)/80);

        vec2 vel = vec2(-pos.x/radius *magn_vel, pos.y/radius* magn_vel) + v;

        return Body({
            avr_mass, 
            pos + center,
            vel
        });
    };

    _s -> loading_queue.add_cluster(cluster);
}

void cosmic_web(Simulation* _s, int nums, double avr_mass){
    noise::module::Perlin perlin; 
    int tot_nums = nums;

    _s -> bodies.reserve(_s -> bodies.size() + nums);

    Cluster2d cluster;
    cluster.number = nums;
    cluster.avr_mass = avr_mass;
    cluster.size = _s -> get_center() *2;
    cluster.name = "cosmic web";

    double g_pull = nums *  _s -> data.G  * avr_mass;

    cluster.builder = [perlin, _s, avr_mass](Cluster2d cl2d, size_t _) {
        vec2 pos = rng_vec2(_s);

        while ((perlin.GetValue(pos.x, pos.y, 0) + 1)/2 < rng_percent()) 
            pos = rng_vec2(_s);

        return Body({
            avr_mass,
            pos,
            vec2({0,0}),
            (perlin.GetValue(pos.x, pos.y, 0) + 1) * 5
        });
    };

    _s -> loading_queue.add_cluster(cluster);
}

/*                   //===---------------------------------------------------------===//
.                    // premade clusters 2d                                           //
.                    //===---------------------------------------------------------===*/


void add_disk3d(Simulation3d* _s, size_t nums, vec3 center, double radius, double thickness, vec3 rotation, double avr_mass = 10e6, vec3 v = {0,0,0}){
    _s -> bodies.reserve(_s -> bodies.size() + nums);

    Cluster3d cluster;
    cluster.number = nums;
    cluster.name = "Disk";

    double g_pull = nums *  _s -> data.G  * avr_mass;

    cluster.builder = [g_pull, radius, v, avr_mass, center, thickness, rotation ](Cluster3d cl3d, size_t _) {
        vec3 pos = rng_point_in_cylinder(radius, 1, thickness);

        double r = pos.sqr_magn();
        double magn_vel = g_pull*r/10000 + std::exp(-(r*r)/80);

        vec3 vel = vec3(-pos.y/r *magn_vel, pos.x/r* magn_vel, rng_percent() - 1) + v;

        pos = rotate_point(pos, rotation.x, rotation.z);
        vel = rotate_point(vel, rotation.x, rotation.z);

        return Body3d({
            avr_mass, 
            pos + center,
            vel
        });
    };

    _s -> loading_queue.add_cluster(cluster);
}


void add_body(Simulation* _s, double mass, vec2 pos, vec2 vel, bool still = false){
    _s -> bodies.push_back(Body(mass, pos, vel));
    _s -> obj ++;
}

void add_body(Simulation* _s, Body b){
    _s -> bodies.push_back(b);
    _s -> obj ++;
}



void add_body(Simulation3d* _s, double mass, vec3 pos, vec3 vel, double theta, double phi){
    _s -> bodies.push_back(Body3d(mass, pos, vel));
    _s -> obj ++;
}

void add_body(Simulation3d* _s, Body3d b){
    _s -> bodies.push_back(b);
    _s -> obj ++;
}


void cosmic_web3d(Simulation3d* _s, int nums, double avr_mass){
    noise::module::Perlin perlin; 
    int tot_nums = nums;

    _s -> bodies.reserve(_s -> bodies.size() + nums);

    Cluster3d cluster;
    cluster.number = nums;
    cluster.avr_mass = avr_mass;
    cluster.size = _s -> get_center() *2;
    cluster.name = "cosmic web";

    double g_pull = nums *  _s -> data.G  * avr_mass;

    cluster.builder = [perlin, _s, avr_mass](Cluster3d cl3d, size_t _) {
        vec3 pos = rng_vec3(_s);

        while ((perlin.GetValue(pos.x, pos.y, pos.z) + 1)/2 < rng_percent()) 
            pos = rng_vec3(_s);

        return Body3d({
            avr_mass,
            pos,
            vec3({0,0,0}),
            (perlin.GetValue(pos.x, pos.y, pos.z) + 1) * 5
        });
    };

    _s -> loading_queue.add_cluster(cluster);
}

}