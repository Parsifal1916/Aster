#pragma once

#include <map>
#include <stdexcept>
#include <string>
#include <chrono>
#include <cassert>
#include <iomanip>

#define CL_TARGET_OPENCL_VERSION 300
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif

#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include "Aster/impl/config.h"
#include "Aster/simulations/basic.h"
#include "Aster/building-api/builder.h"
#include "Aster/graphs/graph_collection.h"

namespace Aster{

template class Simulation<vec2>;
template class Simulation<vec3>;

template <typename T>
bool Simulation<T>::is_fine(){
    bool retval = true;

    for (int i = 0; i < this -> bodies.positions.size(); ++i){
        if (this -> bodies.positions[i].is_fine() && 
            this -> bodies.velocities[i].is_fine() &&
            this -> bodies.accs[i].is_fine())
        continue;

        retval = false;
        break;
    }

    warn_if(!retval, "is_fine() call: simulation is not fine, NaN found in particle's data");

    return retval;
}

/*
* set screen's height, width and bg color (optional)
* @param w_: screen's height
* @param h_: screen's width
* @param c_: screen's bg color
*/
template <typename T>
Simulation<T>* Simulation<T>::set_screen_size(REAL w_, REAL h_, REAL d_){
    data.size = {
        w_, h_, d_
    };

    return this;
}


template <typename T>
struct CoMReducer {
    const std::vector<REAL>& m;                 
    const std::vector<T>& p;

    REAL m_sum {};   
    T mp_sum;

    CoMReducer(const std::vector<REAL>& masses,
               const std::vector<T>& positions)
        : m{masses}, p{positions} {}

    CoMReducer(CoMReducer& other, tbb::split)
        : m{other.m}, p{other.p} {}

    void operator()(const tbb::blocked_range<std::size_t>& r) {
        REAL local_m_sum = m_sum;
        T local_mp_sum = mp_sum;

        for (std::size_t i = r.begin(); i != r.end(); ++i) {
            const REAL mi = static_cast<REAL>(m[i]);
            local_m_sum += mi;
            local_mp_sum += mi * p[i];
        }
        m_sum   = local_m_sum;
        mp_sum  = local_mp_sum;
    }

    void join(const CoMReducer& rhs) {
        m_sum   += rhs.m_sum;
        mp_sum += rhs.mp_sum;
    }
};

template <typename T>
T Simulation<T>::get_center_of_mass(){
    CoMReducer<T> reducer{this -> bodies.masses, this -> bodies.positions};
    tbb::parallel_reduce(tbb::blocked_range<std::size_t>(0, this -> bodies.positions.size()),
                         reducer);

    const double invTotalMass = 1.0 / reducer.m_sum;
    return reducer.mp_sum[0] * invTotalMass;
}


template <typename T>
Simulation<T>* Simulation<T>::set_dt(float dt_){
    assert(dt_ >= 0);
    data.dt = dt_;
    return this;
}

template <typename T>
Simulation<T>* Simulation<T>::set_omega_m(float om_){
    assert(om_ > 0);
    data.omega_m = om_;
    return this;
}

template <typename T>
Simulation<T>* Simulation<T>::set_omega_l(float ol_){
    assert(ol_ > 0);
    data.omega_l = ol_;
    return this;
}

template <typename T>
Simulation<T>* Simulation<T>::set_hubble(float h_){
    assert(h_ > 0);
    data.H_0 = h_;
    return this;
}

template <typename T>
Simulation<T>* Simulation<T>::set_vacuum_density(float d_){
    assert(d_ > 0);
    data.vacuum_density = d_;
    data.root_omega = sqrt(d_);
    return this;
}

template <typename T>
Simulation<T>* Simulation<T>::integrate(size_t times){
    if (times == 0) return this;

    auto start = std::chrono::high_resolution_clock::now();

    for (size_t i = 0; i < times; ++i)
        this -> step();

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<REAL, std::milli> lasted_for = end - start;
    std::ostringstream per_step;
    std::ostringstream total;

    per_step <<  std::fixed << std::setprecision(3) << lasted_for.count() / times;
    total    <<  std::fixed << std::setprecision(3) << lasted_for.count();

    log_info(std::string("Integration report:\n") 
        + std::string("[ + ] Total time:        ") + total.str() + std::string("ms\n")
        + std::string("[ + ] Time for one step: ") + per_step.str()  + std::string("ms\n"));

    return this;
}


template <typename T>
Simulation<T>* Simulation<T>::set_max_frames(unsigned int f_){
    data.max_frames = f_;
    return this;
}

template <typename T>
Simulation<T>* Simulation<T>::set_sim_type(short type){
    assert(type < 4);
    data.type = static_cast<simulation_types>(type);
    return this;
}

template <typename T>
Simulation<T>* Simulation<T>::load(){
    if (has_loaded) return this;
    loading_queue.load(this);
    has_loaded = true;
    this -> calculate_total_mass();
    return this;
}

template <typename T>
Simulation<T>* Simulation<T>::get_force_with(force_func<T> p){
    this -> get_force = p;
    this -> force_used = CUSTOM_F;

    return this;
}

template <typename T>
Simulation<T>* Simulation<T>::use_GPU(){
    this -> GPU_on = true;
    this -> update_forces = get_uf<T>(GPU_UF, this -> force_used);
    return this;
}

template <typename T> 
Simulation<T>* Simulation<T>::update_with(func_ptr<T> p){
    this -> update_bodies = p;
    this -> update_used = CUSTOM_U;
    return this;
}

template < >
inline Simulation<vec3>* Simulation<vec3>::add_graph(typename Graphs::Graph<vec3>::collector_fptr listener, graph_type type){
    assert(type == BETWEEN && "cannot assign this specific function to anything other then a BETWEEN graph");
    this -> between_graphs.push_back({this, listener, type});
    this -> between_graphs.back().name = "Graph" + std::to_string(int(this -> graphs.size()));
    return this;
}

template <>
inline Simulation<vec2>* Simulation<vec2>::add_graph(typename Graphs::Graph<vec2>::collector_fptr listener, graph_type type){
    assert(type == BETWEEN && "cannot assign this specific function to anything other then a BETWEEN graph");
    this -> between_graphs.push_back({this, listener, type});
    this -> between_graphs.back().name = "Graph" + std::to_string(int(this -> graphs.size()));
    return this;
}

template <typename T>
inline Simulation<T>* Simulation<T>::add_graph(typename Graphs::Graph<T>::listener_fptr listener, graph_type type){
    assert(type != BETWEEN && "cannot assign this specific function to be BETWEEN graph");
    this -> graphs.push_back({this, listener, type});
    this -> graphs.back().name = "Graph" + std::to_string(int(this -> graphs.size()));
    return this;
}

template <typename T>
Simulation<T>* Simulation<T>::set_heat_capacity(REAL c_){
    assert(c_ && "heat capacity cannot be 0");
    data.avr_heat_capacity = c_;
    return this;
}

template <typename T>
REAL Simulation<T>::get_height() const{
    return this -> data.size.y * get_scale();
}

template <typename T>
bool Simulation<T>::uses_GPU() const{
    return this -> GPU_on;
}

template <typename T>
REAL Simulation<T>::get_width() const{
    return this -> data.size.x * get_scale();
}

template <typename T>
REAL Simulation<T>::get_depth() const{
    return this -> data.size.z * get_scale();
}


template <typename T>
REAL Simulation<T>::get_render_height() const{
    return this -> data.size.y;
}

template <typename T>
REAL Simulation<T>::get_render_width() const{
    return this -> data.size.x;
}

template <typename T>
REAL Simulation<T>::get_render_depth() const{
    return this -> data.size.z;
}

template <typename T>
int Simulation<T>::get_cores() const{
    return this -> data.NUM_THREADS;
};

template <typename T>
REAL Simulation<T>::get_G() const{
    return this -> data.G;
};

template <typename T>
simulation_types Simulation<T>::get_type() const{
    return this -> data.type;
};


template <typename T>
REAL Simulation<T>::get_c() const{
    return this -> data.c;
}

template <typename T>
REAL Simulation<T>::get_c_sqr() const{
    return this -> data.c_squared;
}

template <typename T>
REAL Simulation<T>::get_dt() const{
    return this -> data.dt;
}

template <typename T>
REAL Simulation<T>::get_e_sqr() const{
    return this -> data.e_squared;
}

template <typename T>
REAL Simulation<T>::get_takeover() const{
    return this -> data.takeover;
}

template <typename T>
REAL Simulation<T>::get_scale() const{
    return this -> data.simulation_scale;
}

template <typename T>
Simulation<T>* Simulation<T>::set_scale(REAL s){
    data.simulation_scale = s;
    return this;
}

template <typename T>
REAL Simulation<T>::get_heat_capacity() const{
    return this -> data.avr_heat_capacity;
}

template <typename T>
REAL Simulation<T>::get_boltzmann() const{
    return this -> data.boltzmann;
}

template <typename T>
void Simulation<T>::trigger_all_graphs(){
    for (auto& graph : graphs)
        graph.trigger();
}

template <typename T>
bool Simulation<T>::has_loaded_yet() const {
    return has_loaded;
}

template <>
inline vec2 Simulation<vec2>::get_center() const{
    return {get_width()/2, get_height()/2};
}

template <>
inline vec3 Simulation<vec3>::get_center() const{
    return {
        get_width()  /2,
        get_height() /2,
        get_depth()  /2
    };
}

template <typename T>
REAL Simulation<T>::get_total_mass() const{
    return total_mass;
}

template <typename T>
REAL Simulation<T>::get_adaptive_coeff() const{
    return data.adaptive_coeff;
}

template <typename T>
Simulation<T>* Simulation<T>::set_adaptive_coeff(REAL s){
    data.adaptive_coeff = s;
    return this;
}


template <typename T>
Simulation<T>* Simulation<T>::calculate_total_mass(){
    total_mass = 0;

    for (const auto& num : this -> bodies.masses)
        total_mass += num;
    
    return this;
}

template <>
inline vec2 Simulation<vec2>::get_corner(int n) const{
    assert(n >= 0 && n < 5);

    return {
        data.size.x * (n % 2),
        data.size.y * (n > 2)
    };
}

template <>
inline vec3 Simulation<vec3>::get_corner(int n) const{
    assert(n >= 0 && n < 9);

    return {
        data.size.x * (n % 2),
        data.size.y * (n==2 || n==3 || n==6 || n==7),
        data.size.z * (n > 3)
    };
}

template <typename T>
REAL Simulation<T>::get_time_passed() const {
    return time_passed;
}

template <typename T>
Simulation<T>* Simulation<T>::collect_hamiltonian(){
    this -> add_graph(Graphs::hamiltonian_collector<T>, ONCE);
    return this;
}

template <typename T>
Simulation<T>* Simulation<T>::collect_error(){
    this -> add_graph(Graphs::error_collector<T>, ONCE);
    return this;
}

template <typename T>
Simulation<T>* Simulation<T>::collect_distance(){
    this -> add_graph(Graphs::distance_collector<T>, FOR_EACH);
    return this;
}

template <typename T>
Simulation<T>* Simulation<T>::get_force_with(force_type t){
    this -> data.selected_force = t;
    this -> get_force = get_force_func<T>(t);
    this -> force_used = t;
    if (uses_GPU()) use_GPU(); // this is cursed af, i know....
    return this;
}

template <typename T>
Simulation<T>* Simulation<T>::update_with(update_type t){
    this -> data.selected_update = t;
    this -> update_bodies = get_update_func<T>(t, uses_GPU());
    this -> update_used = t;
    return this;
}

/**
* @brief single core force update function
*/
template <typename T> 
void single_core_fu(Simulation<T>* _s){
    for (int i = 0; i < _s -> bodies.positions.size(); ++i){
        _s -> bodies.get_acc_of(i).reset();
        _s -> update_pair(i);
    }   
}

/**
* @brief multi core force update function
*/
template <typename T> 
void parallel_fu(Simulation<T>* _s){

    _s -> obj = _s -> bodies.positions.size(); 

    tbb::parallel_for(size_t(0), _s -> bodies.positions.size(), [_s](size_t i){_s -> update_pair(i);});


    for (auto& graph : _s -> between_graphs)
        graph.end_batch();
}

/*                   //===---------------------------------------------------------===//
.                    // SINGLE THREAD IMPLEMENTATION                                  //
.                    //===---------------------------------------------------------===*/

template <typename T>
SingleThread<T>::SingleThread(sim_meta m){
    this -> data = m;
    this -> get_force_with(this -> data.selected_force);
    this -> update_with(this -> data.selected_update);
    this -> data.graph_height *= this -> data.HEIGHT;
    this -> update_forces = single_core_fu;
}

template <typename T>
SingleThread<T>::SingleThread(){
    this -> data = sim_meta();
    this -> data.type = LIGHT;
    this -> get_force_with(this -> data.selected_force);
    this -> update_with(this -> data.selected_update);
    this->update_forces = static_cast<void(*)(Simulation<T>*)>(single_core_fu);

    this -> data.graph_height *= this -> data.size.y;
}

template <typename T>
void SingleThread<T>::step(){
    this -> update_forces(this);
    this -> update_bodies(this);

    for (auto& graph : this -> between_graphs)
        graph.end_batch();

    this -> trigger_all_graphs();
    this -> time_passed++;
}


/*                   //===---------------------------------------------------------===//
.                    // MULTI THREAD IMPLEMENTATION                                   //
.                    //===---------------------------------------------------------===*/

template <typename T>
Parallelized<T>::Parallelized(sim_meta m){
    this -> data = m;
    this -> get_force = get_force_func<T>(this -> data.selected_force);
    this -> update_bodies = get_update_func<T>(this -> data.selected_update, this -> uses_GPU());
    this -> data.graph_height *= this -> data.size.y;
    this -> update_forces = parallel_fu;

    this -> threads.reserve(this -> get_cores()); 
    this -> obj = this -> bodies.size();
}

template <typename T>
Parallelized<T>::Parallelized(){
    this -> data = sim_meta();
    this -> data.type = HEAVY;
    this -> get_force = get_force_func<T>(this -> data.selected_force);
    this -> update_bodies = get_update_func<T>(this -> data.selected_update, this -> uses_GPU());
    this -> data.graph_height *= this -> data.size.y;
    this->update_forces = static_cast<void(*)(Simulation<T>*)>(single_core_fu);


    this -> threads.reserve(this -> get_cores()); 
    this -> obj = this -> bodies.positions.size();
}

template <typename T>
void Parallelized<T>::step(){
    this -> update_forces(this);
    this -> update_bodies(this);

    this -> trigger_all_graphs();
    this -> time_passed++;
}

/*
* sets the maximum number of threads for the given simulation
* if the given number is higher than the number of objs
* it will be topped at the number of objs
*/
template <typename T>
Parallelized<T>* Parallelized<T>::set_max_threads(unsigned int t_){
    this -> get_cores() = t_;
    if (t_ == 0)
        this -> get_cores() = 1;
    if (this -> obj < this -> get_cores())
        this -> get_cores() = this -> obj;

        return this;
}

template <typename T> 
Simulation<T>* Simulation<T>::update_forces_with(func_ptr<T> p){
    this -> update_forces = p;
    this -> force_update_used = CUSTOM_FU;
    return this;
}

template <typename T> 
Simulation<T>* Simulation<T>::update_forces_with(forces_update_type p){
    this -> update_forces = get_uf<T>(p, force_used);
    this -> force_update_used = p;
    return this;
}

template <typename T> 
void update_bundle(Simulation<T>* _s, size_t index){
    unsigned int mult, start, stop;
    mult = _s -> obj/_s -> get_cores();
    start = index * mult;
    stop = (index + 1) * mult;
          
    stop = (stop + mult > _s -> obj) ? _s -> obj : stop;

    for (int i = start; i < stop; ++i)
       _s -> update_pair(i);
}


}


