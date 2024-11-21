#pragma once 

#include <SFML/Graphics.hpp>
#include "../physics/body.h"
#include "simulation.h"

namespace Simulation{
namespace presets{

Body make_spawn_box(unsigned int index){
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
        vel,
        index,
        sf::Color::White//r < radius / 4 ? sf::Color::White : sf::Color::Black
    );
    p.temp = 10;
    if (std::rand() % 2){ p = Body(
        mass, // mass
        pos + vec2({WIDTH-radius, HEIGHT-radius}),
        vel,
        index,
        sf::Color::White
    );}
    
    
    return p; 
}

} //fine presets

    void populate(){
        for (int i = 0; i < obj; ++i) 
            bodies.push_back(presets::make_spawn_box(i));
    }
}//fine sim
