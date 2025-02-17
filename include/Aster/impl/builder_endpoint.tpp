#pragma once

#include <map>
#include <stdexcept>
#include <string>
#include <cassert>

#include "Aster/simulations/basic.h"

#include "Aster/building-api/builder.h"

#include "Aster/graphs/graph_collection.h"

namespace Aster{

template class Simulation<vec2>;
template class Simulation<vec3>;


/*
* set screen's height, width and bg color (optional)
* @param w_: screen's height
* @param h_: screen's width
* @param c_: screen's bg color
*/
template <typename T>
Simulation<T>* Simulation<T>::set_screen_size(double w_, double h_, double d_){
    data.size = {
        w_, h_, d_
    };

    return this;
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
    assert(!has_loaded && "the simulation has already been loaded");
    loading_queue.load(this);
    has_loaded = true;
    return this;
}

template <typename T>
Simulation<T>* Simulation<T>::get_force_with(force_func<T> p){
    this -> get_force = p;
    return this;
}

template <typename T> 
Simulation<T>* Simulation<T>::update_with(func_ptr<T> p){
    this -> update_body = p;
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
Simulation<T>* Simulation<T>::set_heat_capacity(double c_){
    assert(c_ && "heat capacity cannot be 0");
    data.avr_heat_capacity = c_;
    return this;
}

template <typename T>
double Simulation<T>::get_height() const{
    return this -> data.size.y * get_scale();
}

template <typename T>
double Simulation<T>::get_width() const{
    return this -> data.size.x * get_scale();
}

template <typename T>
double Simulation<T>::get_depth() const{
    return this -> data.size.z * get_scale();
}

template <typename T>
double Simulation<T>::get_render_height() const{
    return this -> data.size.y;
}

template <typename T>
double Simulation<T>::get_render_width() const{
    return this -> data.size.x;
}

template <typename T>
double Simulation<T>::get_render_depth() const{
    return this -> data.size.z;
}


template <typename T>
int Simulation<T>::get_cores() const{
    return this -> data.NUM_THREADS;
};

template <typename T>
double Simulation<T>::get_G() const{
    return this -> data.G;
};

template <typename T>
simulation_types Simulation<T>::get_type() const{
    return this -> data.type;
};


template <typename T>
double Simulation<T>::get_c() const{
    return this -> data.c;
}

template <typename T>
double Simulation<T>::get_c_sqr() const{
    return this -> data.c_squared;
}

template <typename T>
double Simulation<T>::get_dt() const{
    return this -> data.dt;
}

template <typename T>
double Simulation<T>::get_e_sqr() const{
    return this -> data.e_squared;
}

template <typename T>
double Simulation<T>::get_takeover() const{
    return this -> data.takeover;
}

template <typename T>
double Simulation<T>::get_scale() const{
    return this -> data.simulation_scale;
}

template <typename T>
Simulation<T>* Simulation<T>::set_scale(double s){
    data.simulation_scale = s;
    return this;
}

template <typename T>
double Simulation<T>::get_heat_capacity() const{
    return this -> data.avr_heat_capacity;
}

template <typename T>
double Simulation<T>::get_boltzmann() const{
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
double Simulation<T>::get_time_passed() const {
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
Simulation<T>* Simulation<T>::get_force_with(force_type t){
    this -> data.selected_force = t;
    this -> get_force = get_force_func<T>(t);
    return this;
}

template <typename T>
Simulation<T>* Simulation<T>::update_with(update_type t){
    this -> data.selected_update = t;
    this -> update_body = get_update_func<T>(t);
    return this;
}

/*                   //===---------------------------------------------------------===//
.                    // SINGLE THREAD IMPLEMENTATION                                  //
.                    //===---------------------------------------------------------===*/

template <typename T>
SingleThread<T>::SingleThread(sim_meta m){
    this -> data = m;
    this -> get_force = get_force_func<T>(this -> data.selected_force);
    this -> update_body = get_update_func<T>(this -> data.selected_update);
    this -> data.graph_height *= this -> data.HEIGHT;
}

template <typename T>
SingleThread<T>::SingleThread(){
    this -> data = sim_meta();
    this -> data.type = LIGHT;
    this -> get_force = get_force_func<T>(this -> data.selected_force);
    this -> update_body = get_update_func<T>(this -> data.selected_update);
    this -> data.graph_height *= this -> data.size.y;
}

template <typename T>
void SingleThread<T>::step(){
    T prev;
    for (Body<T>& body : this -> bodies){  
        prev = body.acceleration;
        body.acceleration.reset();
        this -> update_pair(&body);
        this -> update_body(&body, this );
    }

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
    this -> update_body = get_update_func<T>(this -> data.selected_update);
    this -> data.graph_height *= this -> data.size.y;

    this -> threads.reserve(this -> get_cores()); 
    this -> obj = this -> bodies.size();
}

template <typename T>
Parallelized<T>::Parallelized(){
    this -> data = sim_meta();
    this -> data.type = HEAVY;
    this -> get_force = get_force_func<T>(this -> data.selected_force);
    this -> update_body = get_update_func<T>(this -> data.selected_update);
    this -> data.graph_height *= this -> data.size.y;

    this -> threads.reserve(this -> get_cores()); 
    this -> obj = this -> bodies.size();
}

template <typename T>
void Parallelized<T>::step(){
    this -> obj = this -> bodies.size(); 
    for (int i = 0; i < this -> get_cores(); ++i)
        this -> threads.push_back(std::thread(update_bundle<T>, this, i));

    for (auto& t : threads)
        t.join();

    this -> threads.clear();

    for (auto& graph : this -> between_graphs)
        graph.end_batch();

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
void update_bundle(Simulation<T>* _s, unsigned short index){
    unsigned int mult, start, stop;
    mult = _s -> obj/_s -> get_cores();
    start = index * mult;
    stop = (index + 1) * mult;
          
    stop = (stop + mult > _s -> obj) ? _s -> obj : stop;

    for (int i = start; i < stop; ++i){  
       Body<T>* body = &_s -> bodies[i];
       body -> acceleration.reset();
       _s -> update_pair(body);
       _s -> update_body(body, _s);
    }
}




}


