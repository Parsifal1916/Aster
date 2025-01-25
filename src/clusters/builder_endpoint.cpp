#include <map>
#include <stdexcept>

#include "Aster/building-api/3d_builder.h"
#include "Aster/building-api/builder.h"

#include "Aster/simulations/basic.h"
#include "Aster/simulations/3d_sim_obj.h"
#include "Aster/simulations/sim_obj.h"
#include "Aster/simulations/barnes-hut.h"
#include "Aster/simulations/barnes-hut3d.h"
#include "Aster/simulations/BHT_sim.h"
#include "Aster/simulations/sim_helper.h"

#include "Aster/physics/tool-chain.h"

#include "Aster/graphs/graph_collection.h"

namespace Aster{

//===---------------------------------------------------------===//
// 3d building                                                   //
//===---------------------------------------------------------===//

std::map<std::string, func_ptr3d> update_funcs_3d = {
    {"Euler", update_euler_3d} , 
    {"Leapfrog", update_leapfrog_3d},
    {"Symplectic4", update_symplectic4_3d}
};

std::map<std::string, force_func3d> force_funcs_3d = {
    {"Newton", newtonian_3d} , 
    {"Pn2", pn2_3d}
};

/*
* set screen's height, width and bg color (optional)
* @param w_: screen's height
* @param h_: screen's width
* @param c_: screen's bg color
*/
Simulation3d* Simulation3d::set_screen_size(unsigned int w_, unsigned int h_){
    data.HEIGHT = h_;
    data.WIDTH = w_;
    return this;
}

bool Simulation3d::has_loaded_yet() const {
    return has_loaded;
}

double Simulation3d::get_time_passed() const {
    return time_passed;
}

Simulation3d* Simulation3d::set_dt(float dt_){
    assert(dt_ > 0);
    data.dt = dt_;
    return this;
}

Simulation3d* Simulation3d::set_omega_m(float om_){
    assert(om_ > 0);
    data.omega_m = om_;
    return this;
}

Simulation3d* Simulation3d::set_omega_l(float ol_){
    assert(ol_ > 0);
    data.omega_l = ol_;
    return this;
}

Simulation3d* Simulation3d::set_hubble(float h_){
    assert(h_ > 0);
    data.H_0 = h_;
    return this;
}

Simulation3d* Simulation3d::set_vacuum_density(float d_){
    assert(d_ > 0);
    data.vacuum_density = d_;
    data.root_omega = sqrt(d_);
    return this;
}

Simulation3d* Simulation3d::set_max_frames(unsigned int f_){
    data.max_frames = f_;
    return this;
}

Simulation3d* Simulation3d::set_sim_type(short type){
    assert(type < 4);
    data.type = static_cast<simulation_types>(type);
    return this;
}

Simulation3d* Simulation3d::load(){
    assert(!has_loaded && "the simulation has already been loaded");
    loading_queue.load(this);
    has_loaded = true;
    return this;
}

Simulation3d* Simulation3d::add_graph(Graphs::Graph3d::listener3d_fptr listener, bool for_each_body){
    this -> graphs.push_back({this, listener, for_each_body});
    return this;
}

void Simulation3d::trigger_all_graphs(){
    for (auto& graph : graphs)
        graph.trigger();
}




SingleThread3d::SingleThread3d(sim3d_meta m){
    this -> data = m;
    get_force = force_funcs_3d[data.selected_force];
    update_body = update_funcs_3d[data.selected_update];
    data.graph_height *= data.HEIGHT;
}

SingleThread3d::SingleThread3d(){
    this -> data = sim3d_meta();
    get_force = force_funcs_3d[data.selected_force];
    update_body = update_funcs_3d[data.selected_update];
    data.graph_height *= data.HEIGHT;
}

void SingleThread3d::step(){
    for (auto& body : this -> bodies){  
        update_pair(&body);
        update_body(&body, this );

        if (body.temp >  this -> max_temp) this -> max_temp = body.temp;
        body.temp /= this -> max_temp;
    }

    if (this -> lagrangian > this -> highest_lagrangian) this -> highest_lagrangian = this -> lagrangian;
        
    this -> lagrangians.push_back(this -> lagrangian);
        
    if (this -> lagrangians.size() > data.WIDTH)
        this -> lagrangians.erase(this -> lagrangians.begin());
        
    data.max_temp = 1;
    this -> lagrangian = 0;

    this -> trigger_all_graphs();
    time_passed++;
}

Parallelized3d::Parallelized3d(sim3d_meta m){
    this -> data = m;
    get_force = force_funcs_3d[data.selected_force];
    update_body = update_funcs_3d[data.selected_update];
    data.graph_height *= data.HEIGHT;

    threads.reserve(this -> data.NUM_THREADS); 
    obj = bodies.size();
}

Parallelized3d::Parallelized3d(){
    this -> data = sim3d_meta();
    get_force = force_funcs_3d[data.selected_force];
    update_body = update_funcs_3d[data.selected_update];
    data.graph_height *= data.HEIGHT;

    threads.reserve(this -> data.NUM_THREADS); 
    obj = bodies.size();
}

void Parallelized3d::step(){
    this -> obj = bodies.size(); 
    for (int i = 0; i < this -> data.NUM_THREADS; ++i)
        this -> threads.push_back(std::thread(update_bundle, this, i));

    for (auto& t : threads)
        t.join();

    this -> threads.clear();

    this -> trigger_all_graphs();
    time_passed++;
}
    
/*
* sets the maximum number of threads for the given simulation
* if the given number is higher than the number of objs
* it will be topped at the number of objs
*/
Parallelized3d* Parallelized3d::set_max_threads(unsigned int t_){
    data.NUM_THREADS = t_;
    if (t_ == 0)
        data.NUM_THREADS = 1;
    if (obj < data.NUM_THREADS)
        data.NUM_THREADS = obj;

    return this;
}

void Parallelized3d::update_bundle(Simulation3d* _s, unsigned short index){
    unsigned int mult, start, stop;
    mult = _s -> obj/_s -> data.NUM_THREADS;
    start = index * mult;
    stop = (index + 1) * mult;

    stop = (stop + mult > _s -> obj) ? _s -> obj : stop;
          
    for (int i = start; i < stop; ++i){  
       Body3d* body = &_s -> bodies[i];
       _s -> update_pair(body);
       _s -> update_body(body, _s);

       if (body -> temp > _s -> max_temp) _s -> max_temp = body -> temp;
       body -> temp /= _s -> max_temp;
    }
}

vec3 Simulation3d::get_center() const {
    return {
        data.WIDTH/2,
        data.HEIGHT/2,
        data.depth/2
    };
}

vec3 Simulation3d::get_corner(int n) const {
    return {
        data.WIDTH * (n % 2),
        data.HEIGHT *(n < 2 || (n > 3 && n < 5)),
        data.depth * (n < 3)
    };
}


//===---------------------------------------------------------===//
// 2d building                                                   //
//===---------------------------------------------------------===//


std::uniform_real_distribution<double> angle_rnd(0.0f, 360.0f);
std::uniform_real_distribution<double> normalized_rnd(0.0f, 1.f);

std::map<std::string, func_ptr> update_funcs = {
    {"Euler", update_euler} , 
    {"Leapfrog", update_leapfrog},
    {"Symplectic4", update_symplectic4}
};

std::map<std::string, force_func> force_funcs = {
    {"Newton", newtonian} , 
    {"Pn2", pn2}
};


/*
* set screen's height, width and bg color (optional)
* @param w_: screen's height
* @param h_: screen's width
* @param c_: screen's bg color
*/
Simulation* Simulation::set_screen_size(unsigned int w_, unsigned int h_){
    data.HEIGHT = h_;
    data.WIDTH = w_;
    return this;
}

Simulation* Simulation::set_dt(float dt_){
    assert(dt_ >= 0);
    data.dt = dt_;
    return this;
}

Simulation* Simulation::set_omega_m(float om_){
    assert(om_ > 0);
    data.omega_m = om_;
    return this;
}

Simulation* Simulation::set_omega_l(float ol_){
    assert(ol_ > 0);
    data.omega_l = ol_;
    return this;
}

Simulation* Simulation::set_hubble(float h_){
    assert(h_ > 0);
    data.H_0 = h_;
    return this;
}

Simulation* Simulation::set_vacuum_density(float d_){
    assert(d_ > 0);
    data.vacuum_density = d_;
    data.root_omega = sqrt(d_);
    return this;
}

Simulation* Simulation::set_max_frames(unsigned int f_){
    data.max_frames = f_;
    return this;
}

Simulation* Simulation::set_sim_type(short type){
    assert(type < 4);
    data.type = static_cast<simulation_types>(type);
    return this;
}

Simulation* Simulation::load(){
    assert(!has_loaded && "the simulation has already been loaded");
    loading_queue.load(this);
    has_loaded = true;
    return this;
}

Simulation* Simulation::add_graph(Graphs::Graph2d::listener2d_fptr listener, bool for_each_body){
    this -> graphs.push_back({this, listener, for_each_body});
    return this;
}

void Simulation::trigger_all_graphs(){
    for (auto& graph : graphs)
        graph.trigger();
}



bool Simulation::has_loaded_yet() const {
    return has_loaded;
}


vec2 Simulation::get_center() const{
    return {
        data.WIDTH/2,
        data.HEIGHT/2
    };
}

vec2 Simulation::get_corner(int n) const{
    assert(n >= 0 && n < 5);

    return {
        data.WIDTH * (n % 2),
        data.HEIGHT * (n > 2)
    };
}

double Simulation::get_time_passed() const {
    return time_passed;
}

SingleThread::SingleThread(sim_meta m){
    this -> data = m;
    get_force = force_funcs[data.selected_force];
    update_body = update_funcs[data.selected_update];
    data.graph_height *= data.HEIGHT;
}

SingleThread::SingleThread(){
    this -> data = sim_meta();
    get_force = force_funcs[data.selected_force];
    update_body = update_funcs[data.selected_update];
    data.graph_height *= data.HEIGHT;
}

void SingleThread::step(){
    vec2 prev = {0,0};
    for (auto& body : this -> bodies){  
        prev = body.acceleration;
        update_pair(&body);
        update_body(&body, this );

        if (body.temp >  this -> max_temp) this -> max_temp = body.temp;
        body.temp /= this -> max_temp;
    }

    if (this -> lagrangian > this -> highest_lagrangian) this -> highest_lagrangian = this -> lagrangian;
    
    this -> lagrangians.push_back(this -> lagrangian);
    
    if (this -> lagrangians.size() > data.WIDTH)
        this -> lagrangians.erase(this -> lagrangians.begin());
    
    data.max_temp = 1;
    this -> lagrangian = 0;

    this -> trigger_all_graphs();
    time_passed++;
}

Parallelized::Parallelized(sim_meta m){
    this -> data = m;
    get_force = force_funcs[data.selected_force];
    update_body = update_funcs[data.selected_update];
    data.graph_height *= data.HEIGHT;

    threads.reserve(this -> data.NUM_THREADS); 
    obj = bodies.size();
}

Parallelized::Parallelized(){
    this -> data = sim_meta();
    get_force = force_funcs[data.selected_force];
    update_body = update_funcs[data.selected_update];
    data.graph_height *= data.HEIGHT;

    threads.reserve(this -> data.NUM_THREADS); 
    obj = bodies.size();
}

void update_bundle(Simulation* _s, unsigned short index){
    unsigned int mult, start, stop;
    mult = _s -> obj/_s -> data.NUM_THREADS;
    start = index * mult;
    stop = (index + 1) * mult;
          
    stop = (stop + mult > _s -> obj) ? _s -> obj : stop;

    for (int i = start; i < stop; ++i){  
       Body* body = &_s -> bodies[i];
       _s -> update_pair(body);
       _s -> update_body(body, _s);
       if (body -> temp > _s -> max_temp) _s -> max_temp = body -> temp;
       body -> temp /= _s -> max_temp;
    }
}

void Parallelized::step(){
    this -> obj = bodies.size(); 
    for (int i = 0; i < this -> data.NUM_THREADS; ++i)
        this -> threads.push_back(std::thread(update_bundle, this, i));

    for (auto& t : threads)
        t.join();

    this -> threads.clear();

    this -> trigger_all_graphs();
    time_passed++;
}
    
/*
* sets the maximum number of threads for the given simulation
* if the given number is higher than the number of objs
* it will be topped at the number of objs
*/
Parallelized* Parallelized::set_max_threads(unsigned int t_){
    data.NUM_THREADS = t_;
    if (t_ == 0)
        data.NUM_THREADS = 1;
    if (obj < data.NUM_THREADS)
        data.NUM_THREADS = obj;

        return this;
}



//===---------------------------------------------------------===//
// baking                                                        //
//===---------------------------------------------------------===//


Simulation* bake(simulation_types s){
    switch (s){
     
    case LIGHT:
        return new SingleThread();
     
    case HEAVY:
        return new Parallelized(); 

    case BARNES_HUT:
        return new Barnes::Barnes_Hut();

    case BH_termal:
        return new Barnes::BHT();

    default:
        throw std::runtime_error("Invalid Simulation Type");
        return nullptr;
    }
}

Simulation3d* bake3d(simulation_types s){
    switch (s){
    
    case LIGHT:
        return new SingleThread3d();
    
    case HEAVY:
        return new Parallelized3d(); 

    case BARNES_HUT:
        return new Barnes::Barnes_Hut3d();   

    default:
        break;  
    }

    throw std::runtime_error("Invalid Simulation Type");
    return nullptr;
}

}


