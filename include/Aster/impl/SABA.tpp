#pragma once

#include <cmath>

#include "Aster/physics/tool-chain.h"
#include "Aster/physics/vectors.h"

#include "Aster/simulations/sim_obj.h"
#include "Aster/simulations/basic.h"

#include "Aster/impl/config.h"

// NOTE: i don't trust loop unrolling

namespace Aster{ 

template <typename T> 
void update_SABA1(Simulation<T>* _s){
    update_euler(_s);
}

template <typename T> 
void update_SABA2(Simulation<T>* _s){
    constexpr double c1 = .5 - std::sqrt(3) / 6;
    constexpr double c2 = std::sqrt(3) / 3;
    constexpr double d1 = .5;
    const double dt = _s -> get_dt();

    for_each_body(_s, [c1, dt](Body<T>& body){
        body.position += body.velocity * c1 * dt;
    });

    _s -> update_forces();

    for_each_body(_s, [d1, c2, dt](Body<T>& body){
        body.velocity += body.acceleration * d1 * dt;
        body.position += body.velocity * c2 * dt;        
    });

    _s -> update_forces();

    for_each_body(_s, [d1, c1, dt](Body<T>& body){
        body.velocity += body.acceleration * d1 * dt;
        body.position += body.velocity * c1 * dt;
    });
}

template <typename T> 
void update_SABA3(Simulation<T>* _s){
    constexpr double c1 = .5 - std::sqrt(15) / 10;
    constexpr double c2 = std::sqrt(15)/10;
    constexpr double d1 = 5.0/18.0;
    constexpr double d2 = 4.0/9.0;

    const double dt = _s -> get_dt();

    for_each_body(_s, [c1, dt](Body<T>& body){
        body.position += body.velocity * c1 * dt;
    });

    _s -> update_forces();

    for_each_body(_s, [d1, c2, dt](Body<T>& body){
        body.velocity += body.acceleration * d1 * dt;
        body.position += body.velocity * c2 * dt;        
    });

    _s -> update_forces();

    for_each_body(_s, [d2, c2, dt](Body<T>& body){
        body.velocity += body.acceleration * d2 * dt;
        body.position += body.velocity * c2 * dt;
    });

    _s -> update_forces();

    for_each_body(_s, [d1, c1, dt](Body<T>& body){
        body.velocity += body.acceleration * d1 * dt;
        body.position += body.velocity * c1 * dt;
    });
}

template <typename T> 
void update_SABA4(Simulation<T>* _s){
    constexpr double c1 = .5 - std::sqrt(525 + 70 * std::sqrt(30)) / 70.0;
    constexpr double c2 = (std::sqrt(525 + 70 * std::sqrt(30)) - std::sqrt(525 - 70 * std::sqrt(30))) / 70.0;
    constexpr double c3 = std::sqrt(525 - 70 * std::sqrt(30)) / 35.0;

    constexpr double d1 = .25 - std::sqrt(30) / 72.0;
    constexpr double d2 = .25 + std::sqrt(30) / 72.0;

    const double dt = _s -> get_dt();

    for_each_body(_s, [c1, dt](Body<T>& body){
        body.position += body.velocity * c1 * dt;
    });

    _s -> update_forces();

    for_each_body(_s, [d1, c2, dt](Body<T>& body){
        body.velocity += body.acceleration * d1 * dt;
        body.position += body.velocity * c2 * dt;        
    });

    _s -> update_forces();

    for_each_body(_s, [d2, c3, dt](Body<T>& body){
        body.velocity += body.acceleration * d2 * dt;
        body.position += body.velocity * c3 * dt;        
    });

    _s -> update_forces();

    for_each_body(_s, [d2, c2, dt](Body<T>& body){
        body.velocity += body.acceleration * d2 * dt;
        body.position += body.velocity * c2 * dt;        
    });

    _s -> update_forces();

    for_each_body(_s, [d1, c1, dt](Body<T>& body){
        body.velocity += body.acceleration * d1 * dt;
        body.position += body.velocity * c1 * dt;        
    });
}

template <typename T>
void update_SABA5(Simulation<T>* _s) {
    constexpr double c1 = 0.5 - (std::sqrt(490 + 42 * std::sqrt(105)) + std::sqrt(490 - 42 * std::sqrt(105))) / 84.0;
    constexpr double c2 = std::sqrt(490 - 42 * std::sqrt(105)) / 42.0;
    constexpr double c3 = (std::sqrt(490 + 42 * std::sqrt(105)) - std::sqrt(490 - 42 * std::sqrt(105))) / 84.0;

    constexpr double d1 = (322 - 13 * std::sqrt(70)) / 1800.0;
    constexpr double d2 = (322 + 13 * std::sqrt(70)) / 1800.0;
    constexpr double d3 = 64.0 / 225.0;

    const double dt = _s->get_dt();

    for_each_body(_s, [c1, dt](Body<T>& body) {
        body.position += body.velocity * c1 * dt;
    });

    _s->update_forces();

    for_each_body(_s, [d1, c2, dt](Body<T>& body) {
        body.velocity += body.acceleration * d1 * dt;
        body.position += body.velocity * c2 * dt;
    });

    _s->update_forces();

    for_each_body(_s, [d2, c3, dt](Body<T>& body) {
        body.velocity += body.acceleration * d2 * dt;
        body.position += body.velocity * c3 * dt;
    });

    _s->update_forces();

    for_each_body(_s, [d3, c3, dt](Body<T>& body) { 
        body.velocity += body.acceleration * d3 * dt;
        body.position += body.velocity * c3 * dt;
    });

    _s->update_forces();

    for_each_body(_s, [d2, c2, dt](Body<T>& body) { 
        body.velocity += body.acceleration * d2 * dt;
        body.position += body.velocity * c2 * dt;
    });

    _s->update_forces();

    for_each_body(_s, [d1, c1, dt](Body<T>& body) { 
        body.velocity += body.acceleration * d1 * dt;
        body.position += body.velocity * c1 * dt;
    });
}

template <typename T>
void update_SABA6(Simulation<T>* _s) {
}

}