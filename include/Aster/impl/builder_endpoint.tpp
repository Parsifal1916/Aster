#pragma once

#include <map>
#include <stdexcept>
#include <string>
#include <cassert>

#include "Aster/simulations/basic.h"

#include "Aster/physics/tool-chain.h"
#include "Aster/building-api/builder.h"

#include "Aster/graphs/graph_collection.h"

namespace Aster{

template class Simulation<vec2>;
template class Simulation<vec3>;

template <typename T>
extern  std::map<std::string, func_ptr<T>> update_funcs;

template <typename T>
extern std::map<std::string, force_func<T>> force_funcs;

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
Simulation<T>* Simulation<T>::add_graph(typename Graphs::Graph<T>::listener_fptr listener, bool for_each_body){
    this -> graphs.push_back({this, listener, for_each_body});
    this -> graphs.back().name = "Graph" + std::to_string(int(this -> graphs.size()));
    return this;
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
vec2 Simulation<vec2>::get_center() const{
    return {data.size.x, data.size.y};
}

template <>
vec3 Simulation<vec3>::get_center() const{
    return data.size;
}

template <>
vec2 Simulation<vec2>::get_corner(int n) const{
    assert(n >= 0 && n < 5);

    return {
        data.size.x * (n % 2),
        data.size.y * (n > 2)
    };
}

template <>
vec3 Simulation<vec3>::get_corner(int n) const{
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

/*                   //===---------------------------------------------------------===//
.                    // SINGLE THREAD IMPLEMENTATION                                  //
.                    //===---------------------------------------------------------===*/

template <typename T>
SingleThread<T>::SingleThread(sim_meta m){
    this -> data = m;
    this -> get_force = force_funcs<T>[this -> data.selected_force];
    this -> update_body = update_funcs<T>[this -> data.selected_update];
    this -> data.graph_height *= this -> data.HEIGHT;
}

template <typename T>
SingleThread<T>::SingleThread(){
    this -> data = sim_meta();
    this -> get_force = force_funcs<T>[this -> data.selected_force];
    this -> update_body = update_funcs<T>[this -> data.selected_update];
    this -> data.graph_height *= this -> data.size.y;
}

template <typename T>
void SingleThread<T>::step(){
    T prev;
    for (Body<T>& body : this -> bodies){  
        prev = body.acceleration;
        this -> update_pair(&body);
        this -> update_body(&body, this );

        if (body.temp >  this -> max_temp) this -> max_temp = body.temp;
        body.temp /= this -> max_temp;
    }
    
    this -> data.max_temp = .1;

    this -> trigger_all_graphs();
    this -> time_passed++;
}

/*                   //===---------------------------------------------------------===//
.                    // MULTI THREAD IMPLEMENTATION                                   //
.                    //===---------------------------------------------------------===*/

template <typename T>
Parallelized<T>::Parallelized(sim_meta m){
    this -> data = m;
    this -> get_force = force_funcs<T>[this -> data.selected_force];
    this -> update_body = update_funcs<T>[this -> data.selected_update];
    this -> data.graph_height *= this -> data.size.y;

    this -> threads.reserve(this -> data.NUM_THREADS); 
    this -> obj = this -> bodies.size();
}

template <typename T>
Parallelized<T>::Parallelized(){
    this -> data = sim_meta();
    this -> get_force = force_funcs<T>[this -> data.selected_force];
    this -> update_body = update_funcs<T>[this -> data.selected_update];
    this -> data.graph_height *= this -> data.size.y;

    this -> threads.reserve(this -> data.NUM_THREADS); 
    this -> obj = this -> bodies.size();
}

template <typename T>
void Parallelized<T>::step(){
    this -> obj = this -> bodies.size(); 
    for (int i = 0; i < this -> data.NUM_THREADS; ++i)
        this -> threads.push_back(std::thread(update_bundle<T>, this, i));

    for (auto& t : threads)
        t.join();

    this -> threads.clear();

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
    this -> data.NUM_THREADS = t_;
    if (t_ == 0)
        this -> data.NUM_THREADS = 1;
    if (this -> obj < this -> data.NUM_THREADS)
        this -> data.NUM_THREADS = this -> obj;

        return this;
}


template <typename T> 
void update_bundle(Simulation<T>* _s, unsigned short index){
    unsigned int mult, start, stop;
    mult = _s -> obj/_s -> data.NUM_THREADS;
    start = index * mult;
    stop = (index + 1) * mult;
          
    stop = (stop + mult > _s -> obj) ? _s -> obj : stop;

    for (int i = start; i < stop; ++i){  
       Body<T>* body = &_s -> bodies[i];
       _s -> update_pair(body);
       _s -> update_body(body, _s);
    }
}


}


