#pragma once
#include <vector>
#include <cassert>
#include <thread>

#include "Aster/simulations/sim_obj.h"

#include "Aster/physics/vectors.h"
#include "Aster/physics/body.h"

namespace Aster{   

class Composite2d;

template <typename T>
class CompositeBody{
    public:

    
    CompositeBody(Simulation<T>* _s, size_t segments, double mass, double radius, T initial_position, T initial_velocity = T(0));

    int get_segments() const;
    int get_start() const;
    void update_bodies(Composite2d* _s);

    private:
    size_t segments;
    double radius;
    double total_mass;
    double internal_force;
    size_t init_childs = 0;
    double stiffness = 10; 
    double sidelenght = 0;

    Body<T>& get_next(size_t index) const;
    Body<T>& get_prev(size_t index) const;

    std::vector<Body<T>*> bodies;
};

}

#include "Aster/impl/composite_body.tpp"