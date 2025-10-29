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
#include "Aster/physics/tool-chain.h"

#include "Aster/graphs/graph_collection.h"
#include "Aster/building-api/clusters.h"

namespace Aster{

class Simulation;


vec3 newtonian(REAL m1, REAL m2, vec3 v1, vec3 v2, vec3 p1, vec3 p2, Simulation* _s);


void update_euler(Simulation* _s);


class Solver{
    public:
    Solver(Simulation* _s);
    Solver() {}

    solver_type get_type();
    void set_force(force_type _t);
    void set_force(force_func _t);
    virtual void load() {};
    virtual void compute_forces() {};
    void set_bounds(int lower, int higher = -1);
    Simulation* get_s() {return _s;};
    int get_upper_bound() const;
    int get_lower_bound() const;
    int get_range() const;


    protected:
    
    int lower_int_bound = 0;
    int upper_int_bound = -1;
    force_func get_force;
    force_type _t = NEWTON;
    solver_type type = SINGLE_THREAD;
    Simulation* _s = nullptr;
};


class Updater{
    public:
    Updater(Simulation* _s, int _o, update_type _t);

    update_type get_type();
    virtual void update_bodies();
    int get_order();

    protected:
    int order = 0;  
    func_ptr update;
    update_type type = EULER;
    Simulation* _s = nullptr;
};


struct IntegrationReport{
    REAL total_time = 0;
    REAL time_per_step = 0;
    REAL initial_energy = 0;
    REAL final_energy = 0;
    REAL precision = 0;
    bool tested_precision;
};

class Simulation{
    public:
    float current_a = 0;

    Simulation();
    
    BodyArray bodies;

    Simulation* set_screen_size(REAL  w_, REAL h_, REAL d_ = 1000);
    Simulation* set_dt(float dt_);
    Simulation* set_omega_m(float om_);
    Simulation* set_omega_l(float ol_);
    Simulation* set_hubble(float h_);
    Simulation* set_vacuum_density(float d_);
    Simulation* set_max_frames(unsigned int f_);
    Simulation* set_heat_capacity(REAL c_);
    Simulation* set_scale(REAL s);
    Simulation* set_adaptive_coeff(REAL s);

    bool always_read_pos = true;
    Simulation* load();
    Simulation* load_gpu_buffers();
    Simulation* read_bodies_gpu();
    Simulation* add_graph(typename Graphs::Graph::listener_fptr listener, graph_type type = ONCE);
    Simulation* add_graph(typename Graphs::Graph::collector_fptr listener, graph_type type = BETWEEN);

    Simulation* get_force_with(force_type t);
    Simulation* update_with(update_type t);
    Simulation* update_with(Updater* p);

    Simulation* get_force_with(force_func p);

    Simulation* collect_hamiltonian();
    Simulation* collect_error();
    Simulation* collect_distance();

    virtual IntegrationReport integrate(size_t time, bool precision_test = false);
    Simulation* calculate_total_mass();
    
    REAL get_total_mass() const;
    REAL get_adaptive_coeff() const;
    vec3 get_center_of_mass();
    bool is_fine();

    bool has_loaded_yet() const;
    
    vec3 get_center() const;
    vec3 get_corner(int n) const;

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

    int get_cores() const;

    void step();

    force_type force_used = NEWTON;
    solver_type gravity_solver = PARALLEL;
    update_type integrator = EULER;
    int integrator_order = 0;

    std::vector<Graphs::Graph> between_graphs;
    std::vector<Graphs::Graph> graphs;
    
    ClusterQueue loading_queue;
    std::pair<std::string, REAL> loading_meta = {"", 0};
    Solver* solver = nullptr;
    Updater* updater = nullptr;
    cl_mem positions_cl, accs_cl, velocities_cl, masses_cl;
    int N, selected_N;
    double softening = 1e-14, theta = .2;

    protected:
    REAL total_mass = 0;
    sim_meta data;
 
    REAL
        c_squared,
        time_passed = 0,
        max_temp = 10e5
    ;

    void trigger_all_graphs();
    bool has_loaded = false;
    bool GPU_on = false;
};


func_ptr resolve_gpu_updater(update_type, int ord = 0);
func_ptr bake_update_function(update_type _t, int ord = 0);
Solver* bake_solver(Simulation* _s, solver_type _t);

}