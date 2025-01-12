#pragma once 
#define FNL_IMPL

#include <cmath>
#include "Aster/physics/body.h"
#include "Aster/simulations/simulation.h"
#include "Aster/building-api/FastNoiseLite.h"

namespace Aster{
namespace presets{


vec2 rng_vec2(Simulation* _s){
    return {std::rand() % (int)_s -> data.WIDTH, std::rand() % (int)_s -> data.HEIGHT};
}

vec3 rng_vec3(Simulation3d* _s){
    return {std::rand() % (int)_s -> data.WIDTH, std::rand() % (int)_s -> data.HEIGHT, std::rand() % (int)_s -> data.depth};
}

Body make_spawn_box1(unsigned int index){
    int radius = HEIGHT/10;
    
    // genera l'angolo
    int theta = static_cast<double>(std::rand() % 360);
    // genera il raggio
    double r = (static_cast<double>(std::rand() % radius))/(1+.7*sin(theta)*sin(theta));

    vec2 pos = {std::cos(theta)*r, std::sin(theta)*r};
    double magn_vel = 0*(1/(r*r + .1) + r)/100;
    //if (r <= radius/2)  magn_vel *= 1;
    //if (r <= radius/4)  magn_vel = 1.2 *(1/(r*r*r + .1) + r);
    //if (r <= radius/8)  magn_vel *= 1.5;
    //if (r <= radius/16)  magn_vel *= 1.6;

    vec2 vel = vec2(-std::sin(theta) *magn_vel, std::cos(theta)* magn_vel);
    int mass = normalized_rnd(rng)*10e7;
    Body p(
        mass, // mass
        vec2({radius+cos(theta)*r, radius+sin(theta)*r}),
        vel
    );
    p.temp = 10;
    if (std::rand() % 2){ p = Body(
        mass, // mass
        pos + vec2({WIDTH-radius, HEIGHT-radius}),
        vel
    );}
    
    
    return p; 
}

void print_st_progress(int amount, int total){
    float percent = float(amount) / total *100;
    int equals = percent /2;

    std::cout << "\r[";
    for (int i = 0; i < 50; ++i)
        std::cout << (i < equals ? "=" : " ");

    std::cout << "]" << std::fixed << std::setprecision(1) << percent << "%  " << std::flush;

    if (amount >= total)
        std::cout << "\n\n" << std::endl;
}


double rng_percent(){
    return double(std::rand() % 1000) / 1000.f;
}

void cosmic_web(Simulation* _s, int nums, double avr_mass){
    fnl_state noise = fnlCreateState();
    noise.noise_type = FNL_NOISE_OPENSIMPLEX2;
    int tot_nums = nums;

    std::cout << "\n[ * ] Loading preset \"cosmic web\"... \n";

    _s -> bodies.reserve(_s -> bodies.size() + nums);
    while (nums--){
        vec2 pos = rng_vec2(_s);
        while ((fnlGetNoise2D(&noise, pos.x, pos.y) + 1)/2 < rng_percent()) pos = rng_vec2(_s);
        _s -> bodies.emplace_back(
            avr_mass,
            pos,
            vec2({0,0})
        );
        _s -> bodies.back().temp = (fnlGetNoise2D(&noise, pos.x, pos.y) + 1) * 5;
        print_st_progress(tot_nums - nums, tot_nums);
    }
}

Body make_ring_helper(Simulation* _s, int outer, int inner, vec2 center, double g_pull, double avr_mass = 10e6, vec2 v = {0,0}){//
    double r = (std::rand() % (outer - inner)) + inner;
    double theta = (double)(std::rand() % 360); 

    //// make it into radial coordinates (outer = r^2)
    double magn_vel = g_pull*r/1000 + std::exp(-(r*r)/80);
    vec2 vel = vec2(-std::sin(theta) *magn_vel, std::cos(theta)* magn_vel);

    Body p(
        std::pow(rng_percent(), 2) * avr_mass, 
        vec2({std::cos(theta) * r, std::sin(theta) * r/3}) + center,
        vel + v
    );
    p.temp = 10;
    
    return p; 
}

void add_body(Simulation* _s, double mass, vec2 pos, vec2 vel, bool still = false){
    _s -> bodies.push_back(Body(mass, pos, vel));
    _s -> obj ++;
}

void add_body(Simulation* _s, Body b){
    _s -> bodies.push_back(b);
    _s -> obj ++;
}

void add_disk(Simulation* _s, int outer, int inner, vec2 center, int nums, double avr_mass = 10e6, vec2 v = {0,0}){
    _s -> bodies.reserve(_s -> bodies.size() + nums);

    std::cout << "[ * ] Loading preset \"disk\" at (" << center.x << ", " << center.y << ") ... \n";
    double g_pull = nums *  _s -> data.G  * avr_mass;

    for (int i = 0; i < nums; i++){
        print_st_progress(i, nums);
        add_body(_s, make_ring_helper(_s, outer, inner, center, g_pull, avr_mass, v));
    }
}



void add_body3d(Simulation3d* _s, double mass, vec3 pos, vec3 vel, double theta, double phi){
    _s -> bodies.push_back(Body3d(mass, pos, vel));
    _s -> obj ++;
}

void add_body3d(Simulation3d* _s, Body3d b){
    _s -> bodies.push_back(b);
    _s -> obj ++;
}

vec3 rotate_point(vec3 v, double phi, double theta){
    float x1, z1, y1, z2;
    theta *= M_PI/180 ;
    phi   *= M_PI/180 ; 

    x1 = v.x * std::cos(phi) + v.z * std::sin(phi);
    z1 = - v.x * std::sin(phi) + v.z *std::cos(phi);

    y1 = v.y * std::cos(theta) + z1 *std::sin(theta);
    z2 = -v.y *std::sin(theta) + z1 *std::cos(theta);

    return {x1, y1, z2};
}

Body3d make_ring_helper3d(Simulation3d* _s, int outer, int inner, vec3 center, double g_pull, int thickness,  double phi, double theta, double avr_mass = 10e6, vec3 v = {0,0,0}){//
    double r = (std::rand() % (outer - inner)) + inner;
    double t = (double)(std::rand() % 360); 

    //// make it into radial coordinates (outer = r^2)
    double magn_vel =std::sqrt(g_pull*r/10e4);//g_pull*r/10e6- std::exp(-(r*r)/80);
    vec3 vel = vec3(std::sin(t) *magn_vel, -std::cos(t)* magn_vel, thickness ? (std::rand() % (2*thickness)) - thickness : 0);

    vec3 pos = vec3({std::cos(t) * r, std::sin(t) * r/3, thickness ? (std::rand() % (2*thickness)) - thickness : 0});
    pos = rotate_point(pos, phi, theta);
    vel = rotate_point(vel, phi, theta);

    Body3d p(
        std::pow(rng_percent(), 2) * avr_mass, 
        pos + center,
        vel
    );

    p.temp = 10;
    
    return p; 
}

void add_disk3d(Simulation3d* _s, int outer, int inner, vec3 center, int nums, double thickness, double phi, double theta, double avr_mass = 10e8, vec3 v = {0, 0,0}){
    _s -> bodies.reserve(_s -> bodies.size() + nums);
    

    std::cout << "\n[ * ] Loading preset \"disk\" at (" << center.x << ", " << center.y << ") ... \n";
    double g_pull = nums *  _s -> data.G  * avr_mass;
    

    for (int i = 0; i < nums; i++){
        print_st_progress(i, nums);
        add_body3d(_s, make_ring_helper3d(_s, outer, inner, center, g_pull, thickness, phi, theta, avr_mass, v));
    }
}

vec3 rng_3d(int max_r){
    return {
        (rand() % (max_r*2)) - max_r,
        (rand() % (max_r*2)) - max_r,
        (rand() % (max_r*2)) - max_r
    };
}

Body3d pick_3d_point(Simulation3d* _s, int max_r, vec3 center){
    vec3 retval = rng_3d(max_r);

    while (retval.magnitude() > max_r) retval = rng_3d(max_r);

    return Body3d(10e10, retval +center, {0,0,0}); 
}

void rng_sphere(Simulation3d* _s, int max_r, vec3 center, int nums){
    _s -> bodies.reserve(_s -> bodies.size() + nums);

    std::cout << "[ * ] Loading preset \"sphere\" at (" << center.x << ", " << center.y << ") ... \n";

    for (int i = 0; i < nums; i++){
        print_st_progress(i, nums);
        add_body3d(_s, pick_3d_point(_s, max_r, center));
    }
}

void cosmic_web3d(Simulation3d* _s, int nums, double avr_mass){
    fnl_state noise = fnlCreateState();
    noise.noise_type = FNL_NOISE_OPENSIMPLEX2;
    int tot_nums = nums;

    std::cout << "\n[ * ] Loading preset \"cosmic web3d\"... \n";

    _s -> bodies.reserve(_s -> bodies.size() + nums);
    while (nums--){
        vec3 pos = rng_vec3(_s);
        while ((fnlGetNoise3D(&noise, pos.x, pos.y, pos.z) + 1)/2 < rng_percent()) pos = rng_vec3(_s);
        _s -> bodies.emplace_back(
            avr_mass,
            pos,
            vec3({0,0,0})
        );
        _s -> bodies.back().temp = (fnlGetNoise3D(&noise, pos.x, pos.y, pos.z) + 1) * 5;
        print_st_progress(tot_nums - nums, tot_nums);
    }
}


} //fine presets
}//fine sim
