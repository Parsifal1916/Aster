#pragma once
#include <random>
#include <string>

#include "Aster/physics/vectors.h"
#include "Aster/physics/tool-chain.h"

#include "Aster/simulations/basic.h"
#include "Aster/simulations/sim_obj.h"
#include "Aster/simulations/3d_sim_obj.h"

#include "Aster/building-api/clusters.h"
#include "Aster/building-api/sim_meta.h"

namespace Aster{

class Simulation{
    public:
    float current_a = 0;

    std::vector<double> lagrangians;
    std::vector<Body> bodies;

    Simulation* set_numof_objs(unsigned int n_);
    Simulation* set_screen_size(unsigned int w_, unsigned int h_);
    Simulation* set_dt(float dt_);
    Simulation* set_omega_m(float om_);
    Simulation* set_omega_l(float ol_);
    Simulation* set_hubble(float h_);
    Simulation* set_vacuum_density(float d_);
    Simulation* set_max_frames(unsigned int f_);
    Simulation* set_sim_type(short type);
    Simulation* load();
    bool has_loaded_yet() const;

    friend void update_bundle(Simulation*, short unsigned int);
    
    vec2 get_center() const;
    vec2 get_corner(int n) const;
    
    virtual void step(){assert(1);}

    virtual void update_pair(Body* b1){
        for (const auto& b2 : bodies){
            if (&b2 == b1) continue;
    
            b1 -> acceleration += get_force(
                b1 -> mass, b2.mass,
                b1 -> velocity, b2.velocity,
                b1 -> position, b2.position,
                this
            ) / b1 -> mass;
    
        }
        lagrangian += .5 * b1 -> velocity.sqr_magn() * b1 -> mass - b1 -> acceleration.magnitude() * b1 -> mass;
    }

    func_ptr update_body;
    force_func get_force;
    
    sim_meta data;
    int obj;
    
    double lagrangian, highest_lagrangian;
    struct Queue2d loading_queue;
    std::pair<std::string, double> loading_meta = {"", 0};

    protected:
    double
        c_squared,
        time_passed = 0,
        max_temp = 10e5

    ;

    bool has_loaded = false;

};

}
