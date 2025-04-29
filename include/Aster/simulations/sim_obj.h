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
    Simulation<T>* set_heat_capacity(double c_);
    Simulation<T>* set_scale(double s);
    Simulation<T>* set_adaptive_coeff(double s);
    
    Simulation<T>* load();
    Simulation<T>* add_graph(typename Graphs::Graph<T>::listener_fptr listener, graph_type type = ONCE);
    Simulation<T>* add_graph(typename Graphs::Graph<T>::collector_fptr listener, graph_type type = BETWEEN);

    Simulation<T>* get_force_with(force_type t);
    Simulation<T>* update_with(update_type t);
    Simulation<T>* get_force_with(force_func<T> p);
    Simulation<T>* update_with(func_ptr<T> p);

    Simulation<T>* collect_hamiltonian();
    Simulation<T>* collect_error();
    Simulation<T>* collect_distance();

    Simulation<T>* calculate_total_mass();
    double get_total_mass() const;
    double get_adaptive_coeff() const;

    bool has_loaded_yet() const;
    
    T get_center() const;
    T get_corner(int n) const;

    double get_height() const;
    double get_width() const;
    double get_depth() const;

    double get_render_height() const;
    double get_render_width() const;
    double get_render_depth() const;

    double get_G() const;
    double get_time_passed() const;
    double get_c() const;
    double get_c_sqr() const;
    double get_dt() const;
    double get_e_sqr() const;
    double get_takeover() const;
    double get_scale() const;
    double get_heat_capacity() const;
    double get_boltzmann() const;

    Simulation<T>* use_simd();
    bool is_using_simd() const ;

    simulation_types get_type() const;
    
    int get_cores() const;

    virtual void step() {}

    virtual void update_pair(Body<T>* b1){
        for (auto& b2 : bodies){
            if (&b2 == b1) continue;
    
            b1 -> acceleration += get_force(
                b1 -> mass, b2.mass,
                b1 -> velocity, b2.velocity,
                b1 -> position, b2.position,
                this
            ) / b1 -> mass;

            for (auto& graph : between_graphs)
                graph.trigger_on(b1, &b2);
    
        }
    }

    func_ptr<T> update_bodies;
    force_func<T> get_force;
    func_ptr<T> update_forces;
    int obj;
    std::vector<Graphs::Graph<T>> between_graphs;
    ClusterQueue<T> loading_queue;
    std::pair<std::string, double> loading_meta = {"", 0};


    protected:
    double total_mass = 0;
    sim_meta data;
    std::vector<Graphs::Graph<T>> graphs;


    double
        c_squared,
        time_passed = 0,
        max_temp = 10e5
    ;

    void trigger_all_graphs();
    bool has_loaded = false;
    bool simd_on = false;

};

template <typename T> 
void update_bundle(Simulation<T>*, short unsigned int);
}

#include "Aster/impl/builder_endpoint.tpp"
