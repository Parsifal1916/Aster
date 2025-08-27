#pragma once
#include <random>
#include <string>

#define CL_TARGET_OPENCL_VERSION 300
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif

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
    
    BodyArray<T> bodies;

    Simulation<T>* set_numof_objs(unsigned int n_);
    Simulation<T>* set_screen_size(REAL  w_, REAL h_, REAL d_ = 1000);
    Simulation<T>* set_dt(float dt_);
    Simulation<T>* set_omega_m(float om_);
    Simulation<T>* set_omega_l(float ol_);
    Simulation<T>* set_hubble(float h_);
    Simulation<T>* set_vacuum_density(float d_);
    Simulation<T>* set_max_frames(unsigned int f_);
    Simulation<T>* set_sim_type(short type);
    Simulation<T>* set_heat_capacity(REAL c_);
    Simulation<T>* set_scale(REAL s);
    Simulation<T>* set_adaptive_coeff(REAL s);
    
    virtual Simulation<T>* load();
    Simulation<T>* add_graph(typename Graphs::Graph<T>::listener_fptr listener, graph_type type = ONCE);
    Simulation<T>* add_graph(typename Graphs::Graph<T>::collector_fptr listener, graph_type type = BETWEEN);

    Simulation<T>* get_force_with(force_type t);
    Simulation<T>* update_with(update_type t);
    Simulation<T>* get_force_with(force_func<T> p);
    Simulation<T>* update_with(func_ptr<T> p);
    Simulation<T>* update_forces_with(func_ptr<T> p);
    Simulation<T>* update_forces_with(forces_update_type p);

    Simulation<T>* collect_hamiltonian();
    Simulation<T>* collect_error();
    Simulation<T>* collect_distance();

    Simulation<T>* integrate(size_t time);
    Simulation<T>* calculate_total_mass();
    virtual Simulation<T>* use_GPU();
    
    REAL get_total_mass() const;
    REAL get_adaptive_coeff() const;
    T get_center_of_mass();
    bool is_fine();

    bool has_loaded_yet() const;
    
    T get_center() const;
    T get_corner(int n) const;

    REAL get_height() const;
    REAL get_width() const;
    REAL get_depth() const;

    REAL get_render_height() const;
    REAL get_render_width() const;
    REAL get_render_depth() const;

    REAL get_G() const;
    REAL get_time_passed() const;
    REAL get_c() const;
    REAL get_c_sqr() const;
    REAL get_dt() const;
    REAL get_e_sqr() const;
    REAL get_takeover() const;
    REAL get_scale() const;
    REAL get_heat_capacity() const;
    REAL get_boltzmann() const;

    bool uses_GPU() const;

    simulation_types get_type() const;
    
    int get_cores() const;

    virtual void step() {}

    virtual void update_pair(size_t  b1){
        for (int b2 = 0; b2 < bodies.positions.size(); ++b2){
            if (b2 == b1) continue;
    
            this -> bodies.get_acc_of(b1) += get_force(
                bodies.get_mass_of(b1), bodies.get_mass_of(b2),
                bodies.get_velocity_of(b1), bodies.get_velocity_of(b2),
                bodies.get_position_of(b1), bodies.get_position_of(b2),
                this
            ) / bodies.get_mass_of(b1);

            for (auto& graph : between_graphs)
                graph.trigger_on(b1, b2);
    
        }
    }

    func_ptr<T> update_bodies;
    force_func<T> get_force;
    func_ptr<T> update_forces;
    int obj;
    std::vector<Graphs::Graph<T>> between_graphs;
    ClusterQueue<T> loading_queue;
    std::pair<std::string, REAL> loading_meta = {"", 0};


    protected:
    REAL total_mass = 0;
    sim_meta data;

    force_type force_used = NEWTON;
    forces_update_type force_update_used = PARALLEL;
    update_type update_used = EULER;
    
    std::vector<Graphs::Graph<T>> graphs;
 
    REAL
        c_squared,
        time_passed = 0,
        max_temp = 10e5
    ;

    void trigger_all_graphs();
    bool has_loaded = false;
    bool GPU_on = false;
};

template <typename T> 
void update_bundle(Simulation<T>*, short unsigned int);
}

#include "Aster/impl/builder_endpoint.tpp"
