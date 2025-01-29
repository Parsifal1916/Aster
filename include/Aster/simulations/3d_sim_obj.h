#pragma once
#include <random>
#include <cassert>
#include <string>
#include <vector>

#include "Aster/simulations/basic.h"
#include "Aster/building-api/sim_meta.h"
#include "Aster/building-api/clusters.h"
#include "Aster/graphs/graph_collection.h"

namespace Aster{

class Simulation3d{
    public:
    int obj;
    sim3d_meta data;

    double
        c_squared,
        time_passed = 0,
        max_temp = 10e5
    ;

    float current_a = 0;

    std::vector<Body3d> bodies;

    std::uniform_real_distribution<double> get_rndX;
    std::uniform_real_distribution<double> get_rndY;

    Simulation3d* set_numof_objs(unsigned int n_);
    Simulation3d* set_screen_size(unsigned int w_, unsigned int h_);
    Simulation3d* set_dt(float dt_);
    Simulation3d* set_omega_m(float om_);
    Simulation3d* set_omega_l(float ol_);
    Simulation3d* set_hubble(float h_);
    Simulation3d* set_vacuum_density(float d_);
    Simulation3d* set_max_frames(unsigned int f_);
    Simulation3d* set_sim_type(short type);
    Simulation3d* load();

    Simulation3d* add_graph(Graphs::Graph3d::listener3d_fptr listener, bool for_each_body = false);

    bool has_loaded_yet() const;



    vec3 get_center() const;
    vec3 get_corner(int n) const;

    double get_time_passed() const ;

    virtual void step(){assert(1);}

    virtual void update_pair(Body3d* b1){
        for (const auto& b2 : bodies){
            if (&b2 == b1) continue;
    
            b1 -> acceleration += get_force(
                b1 -> mass, b2.mass,
                b1 -> velocity, b2.velocity,
                b1 -> position, b2.position,
                this
            ) / b1 -> mass;
    
        }
    }

    func_ptr3d update_body;
    force_func3d get_force;

    std::pair<std::string, double> loading_meta = {"", 0};
    
    Queue3d loading_queue;

    protected:
    void trigger_all_graphs();
    std::vector<Graphs::Graph3d> graphs;
    bool has_loaded = false;
};

}
