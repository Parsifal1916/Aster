#pragma once

#include <cassert>
#include "Aster/simulations/composite.h"

namespace Aster{

    
template <> 
CompositeBody<vec2>::CompositeBody(Simulation<vec2>* _s, size_t segments, double mass, double radius, vec2 initial_position, vec2 initial_velocity)
    : total_mass(mass), radius(radius), segments(segments){

    this -> bodies.reserve(segments + 1);
    double mass_per_body = mass / (segments + 1);
    init_childs = _s -> bodies.size();

    internal_force = mass * mass_per_body* _s -> get_G() / (radius * radius);
    sidelenght = 2 * M_PI * radius / segments;

    // pushes the center
    _s -> bodies.push_back(Body<vec2>(
        mass_per_body,
        initial_position,
        initial_velocity
    ));

    // saves the center for itself
    this -> bodies.push_back(&_s -> bodies.back());

    vec2 position = {0,0};
    double angle = 0;

    for (int i = 0; i < segments; i++){
        angle = M_PI * 2 * i / segments;
        
        position.x = std::cos(angle) * radius;
        position.y = std::sin(angle) * radius;

        _s -> bodies.push_back(Body<vec2>(
            mass_per_body,
            position + initial_position,
            initial_velocity
        ));

        this -> bodies.push_back(&_s -> bodies.back());
    }
}

template <typename T>
inline Body<T>& CompositeBody<T>::get_next(size_t index) const{
    assert(index < this -> bodies.size() - 1);

    if (index + 1 > this -> bodies.size() - 1)
        return this -> bodies[0];
    return this -> bodies[index + 1];
}

template <typename T>
inline Body<T>& CompositeBody<T>::get_prev(size_t index) const{
    assert(index < this -> bodies.size() - 1);
    
    if (index - 1 < 0)
        return this -> bodies.back();
    return this -> bodies[index - 1];
}

template <typename T>
int CompositeBody<T>::get_segments() const{
    return this -> segments;
}

template <typename T>
int CompositeBody<T>::get_start() const{
    return this -> init_childs;
}

template <typename T>
void apply_spring(Composite2d* _s, Body<T>& b1, Body<T>& b2, double bound, double stiffness, double base = 0) {
    T arrow = b2.position - b1.position;
    double r = arrow.magnitude();

    arrow = (arrow / r * (r - bound) * stiffness + base ) * _s -> get_dt();

    b1.acceleration +=  arrow / b1.mass;
    b2.acceleration += -arrow / b2.mass;
}

template <typename T> 
void CompositeBody<T>::update_bodies(Composite2d* _s){
    for (int i = get_start() + 1; i < get_segments() + get_start(); i++){
        // between pairs of objects
        apply_spring(_s, _s -> bodies[i], _s -> bodies[i+1], sidelenght, stiffness);

        // betweeen the center and object
        apply_spring(_s, _s -> bodies[i], _s -> bodies[get_start()], radius, stiffness, internal_force);
    }

    apply_spring(_s, _s -> bodies[get_segments() + get_start()], _s -> bodies[get_start() + 1], sidelenght, stiffness);
    apply_spring(_s, _s -> bodies[get_segments() + get_start()], _s -> bodies[get_start()], radius, stiffness, internal_force);    
}

Composite2d* add_composite(Composite2d* _s, size_t segments, double mass, double radius, vec2 initial_position, vec2 initial_velocity){
    _s -> composites.push_back(CompositeBody<vec2>(
        _s,
        segments, 
        mass, 
        radius, 
        initial_position, 
        initial_velocity
    ));
    
    return _s;
}

}