#pragma once 

#include <SFML/Graphics.hpp>
#include "../physics/body.h"
#include "simulation.h"

namespace Simulation{
namespace presets{

Body make_ellipse(unsigned int index, vec2 center){
    int radius = HEIGHT/10;
    
    int theta = static_cast<double>(std::rand() % 360);
    double r = (static_cast<double>(std::rand() % radius))/(1+.7*sin(theta)*sin(theta));

    vec2 pos = {std::cos(theta)*r, std::sin(theta)*r};
    double magn_vel = 0*(1/(r*r + .1) + r)/100;

    vec2 vel = vec2(-std::sin(theta) *magn_vel, std::cos(theta)* magn_vel);
    int mass = normalized_rnd(rng)*10e7;
    Body p(
        mass, // mass
        center + vec2({cos(theta)*r, sin(theta)*r}),
        vel
    );
    p.temp = 10;
    
    return p; 
}

Body make_ring_helper(int inner, int outer, vec2 center){
    double x = (std::rand() % (2*outer)) - outer, y = (std::rand() % (2*outer)) - outer;

    while( x*x + y*y > outer*outer || x*x + y*y < outer*outer/4) { // makes sure it's in the circle
        x = (std::rand() % (2*outer)) - outer;
        y = (std::rand() % (2*outer)) - outer;
    }

    // make it into radial coordinates (outer = r^2)
    outer = x*x + y*y;
    double theta = std::atan(x/(y+0.00001));

    double magn_vel = obj*G/(outer);
    vec2 vel = vec2(-std::sin(theta) *magn_vel, std::cos(theta)* magn_vel);

    Body p(
        1, 
        vec2({x, y}) + center,
        vel
    );
    p.temp = 10;
    
    return p; 
}

void make_ring(int outer, int inner, vec2 center, int nums){
    bodies.reserve(bodies.size() + nums);

    while (nums--){
        bodies.push_back(make_ring_helper(outer, inner, center));
    }
}

void add_body(double mass, vec2 pos, vec2 vel, sf::Color col, bool still = false){
    bodies.push_back(Body(mass, pos, vel, col, still));
}

} //fine presets
}//fine sim
