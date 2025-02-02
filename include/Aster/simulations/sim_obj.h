#pragma once
#include <random>
#include <string>

#include "Aster/physics/vectors.h"

#include "Aster/simulations/basic.h"

#include "Aster/building-api/clusters.h"
#include "Aster/building-api/sim_meta.h"

#include "Aster/graphs/graph_collection.h"

namespace Aster{

template <typename T> struct ClusterQueue;

template <typename T>
class Simulation{
    public:
    float current_a = 0;
    
    std::vector<Body<T>> bodies;

    Simulation<T>* set_numof_objs(unsigned int n_);
    Simulation<T>* set_screen_size(double  w_, double h_, double d_ = 1000);
    Simulation<T>* set_dt(float dt_);
    Simulation<T>* set_omega_m(float om_);
    Simulation<T>* set_omega_l(float ol_);
    Simulation<T>* set_hubble(float h_);
    Simulation<T>* set_vacuum_density(float d_);
    Simulation<T>* set_max_frames(unsigned int f_);
    Simulation<T>* set_sim_type(short type);
    Simulation<T>* load();

    Simulation<T>* add_graph(typename Graphs::Graph<T>::listener_fptr listener, bool for_each_body = false);

<<<<<<< HEAD
=======
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
    double get_time_passed();
>>>>>>> main
    bool has_loaded_yet() const;
    
    T get_center() const;
    T get_corner(int n) const;

    double get_time_passed() const ;
    
    virtual void step() {}

    virtual void update_pair(Body<T>* b1){
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

    func_ptr<T> update_body;
    force_func<T> get_force;
    
    sim_meta data;
    int obj;

    ClusterQueue<T> loading_queue;
    std::pair<std::string, double> loading_meta = {"", 0};


    protected:

    std::vector<Graphs::Graph<T>> graphs;
    
    double
        c_squared,
        time_passed = 0,
        max_temp = 10e5
    ;

<<<<<<< HEAD
    void trigger_all_graphs();
=======
    std::vector<struct Graph2d*> graphs;

>>>>>>> main
    bool has_loaded = false;

};

template <typename T> 
void update_bundle(Simulation<T>*, short unsigned int);
}

#include "Aster/impl/builder_endpoint.tpp"
