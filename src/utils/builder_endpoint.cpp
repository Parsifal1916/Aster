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

#include "Aster/building-api/GPU_endpoint.h"
#include "Aster/simulations/sim_obj.h"
#include "Aster/simulations/basic.h"
#include "Aster/building-api/builder.h"
#include "Aster/graphs/graph_collection.h"

namespace Aster{
using namespace tbb;

//===---------------------------------------------------------===//
//                 SOLVERS IMPLEMENTATION                        //
//===---------------------------------------------------------===//

Solver::Solver(Simulation* _s){
    this -> _s = _s;
}

int Solver::get_upper_bound() const{
    return (upper_int_bound < 0) 
        ? this -> _s -> bodies.positions.size()
        : upper_int_bound; 
}

int Solver::get_lower_bound() const{
    return lower_int_bound; 
}

int Solver::get_range() const{
    return get_upper_bound() - get_lower_bound(); 
}


solver_type Solver::get_type(){
    return this -> type;
}


void Solver::set_force(force_func _t){
    this -> get_force = _t;
    this -> _t = CUSTOM_F;  
}


void Solver::set_bounds(int lower, int upper){
    this -> lower_int_bound = lower;
    this -> upper_int_bound = upper;
}

//===---------------------------------------------------------===//
//                 UPDATER IMPLEMENTATION                        //
//===---------------------------------------------------------===//


update_type Updater::get_type(){
    return this -> type;
}


void Updater::update_bodies(){
    this -> update(this -> _s);
}


int Updater::get_order(){
    return this -> order;
};


Updater::Updater(Simulation* s, int _o, update_type _t)
    : type(_t), _s(s), order(_o){
    
    if (s -> uses_GPU())
        this -> update = resolve_gpu_updater(_t, _o);
    else
        this -> update = bake_update_function(_t, _o);
}

//===---------------------------------------------------------===//
//                 SIM-UTILS IMPLEMENTATION                      //
//===---------------------------------------------------------===//

Simulation* Simulation::load_gpu_buffers(){
    using namespace GPU;
    if (!has_initialized) {
        init_opencl();
        has_initialized = true;
    }

    if (warn_if(!this -> GPU_on, "Cannot load gpu buffers without this -> GPU_on = true")) return this;
    size_t N = this -> bodies.positions.size();
    log_info("Loading internal AOB buffers... ");

    positions_cl = create_vbuffer(this -> bodies.positions);
    velocities_cl = create_vbuffer(this -> bodies.velocities);
    accs_cl = create_vbuffer(this -> bodies.accs);
    masses_cl = create_vbuffer(this -> bodies.masses);

    return this;
}

Simulation::Simulation(){

}

bool Simulation::is_fine(){
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

//===---------------------------------------------------------===//
//                 SIM-SETTERS IMPLEMENTATION                    //
//===---------------------------------------------------------===//

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


/*
* set screen's height, width and bg color (optional)
* @param w_: screen's height
* @param h_: screen's width
* @param c_: screen's bg color
*/

Simulation* Simulation::set_screen_size(REAL w_, REAL h_, REAL d_){
    data.size = {
        w_, h_, d_
    };

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

Simulation* Simulation::set_heat_capacity(REAL c_){
    assert(c_ && "heat capacity cannot be 0");
    data.avr_heat_capacity = c_;
    return this;
}

Simulation* Simulation::set_scale(REAL s){
    data.simulation_scale = s;
    return this;
}


//===---------------------------------------------------------===//
//                 SIM-GETTERS IMPLEMENTATION                    //
//===---------------------------------------------------------===//

struct CoMReducer {
    const std::vector<REAL>& m;                 
    const std::vector<vec3>& p;

    REAL m_sum {};   
    vec3 mp_sum;

    CoMReducer(const std::vector<REAL>& masses,
               const std::vector<vec3>& positions)
        : m{masses}, p{positions} {}

    CoMReducer(CoMReducer& other, tbb::split)
        : m{other.m}, p{other.p} {}

    void operator()(const tbb::blocked_range<std::size_t>& r) {
        REAL local_m_sum = m_sum;
        vec3 local_mp_sum = mp_sum;

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

vec3 Simulation::get_center() const{
    return {
        this -> get_width()/2,
        this -> get_height()/2,
        this -> get_depth()/2
    };
}

vec3 Simulation::get_center_of_mass(){
    CoMReducer reducer{this -> bodies.masses, this -> bodies.positions};
    tbb::parallel_reduce(tbb::blocked_range<std::size_t>(0, this -> bodies.positions.size()),
                         reducer);

    const double invTotalMass = 1.0 / reducer.m_sum;
    return reducer.mp_sum[0] * invTotalMass;
}




void Simulation::step(){
    this -> selected_N = this -> solver -> get_range();
    this -> N = this -> bodies.positions.size();
    this -> updater -> update_bodies();

    this -> trigger_all_graphs();
    if (always_read_pos && uses_GPU())
        read_bodies_gpu();
    this -> time_passed++;
}



IntegrationReport Simulation::integrate(size_t times, bool precision_test){
    using namespace Graphs;
    if (times == 0) return {};
    std::ostringstream aft, bef, prec;

    REAL e_before, e_after;

    if (precision_test){
        if (uses_GPU())
            read_bodies_gpu();
        e_before = get_total_energy(this);
    }

    bool prev = this->always_read_pos; 
    this->always_read_pos = false;

    auto start = std::chrono::high_resolution_clock::now();

    for (size_t i = 0; i < times; ++i)
        this -> step();

    auto end = std::chrono::high_resolution_clock::now();

    this->always_read_pos = prev;
    if (precision_test){
        if (uses_GPU())
            read_bodies_gpu();
        e_after = get_total_energy(this);
    }

    std::chrono::duration<REAL, std::milli> lasted_for = end - start;
    std::ostringstream per_step;
    std::ostringstream total;

    per_step <<  std::fixed << std::setprecision(3) << lasted_for.count() / times;
    total    <<  std::fixed << std::setprecision(3) << lasted_for.count();

    std::string msg = std::string("Integration report:\n") 
        + std::string("[ + ] Total time:        ") + total.str() + std::string("ms\n")
        + std::string("[ + ] Time for one step: ") + per_step.str()  + std::string("ms\n");

    REAL acc = std::abs(e_before - e_after);
    if (!e_before) acc = 0;
    else acc /= std::abs(e_before);

    std::string add = "";
    if (precision_test){
        aft << std::scientific << std::setprecision(3) << e_after;
        bef << std::scientific << std::setprecision(3) << e_before;
        prec  << std::scientific << std::setprecision(3) << acc;
        add = std::string("[ + ] Total Energy Before:  ") + std::string(bef.str()  + std::string("J")) 
          + std::string("\n[ + ] Total Energy After:   ") + std::string(aft.str()) + std::string("J")
          + std::string("\n[ + ] Integrator Precision: ") + std::string(prec.str());
    }

    log_info(msg + add);

    return {
        lasted_for.count(), 
        lasted_for.count() / times, 
        e_before, 
        e_after,  
        acc,
        precision_test
    };
}

Simulation* Simulation::read_bodies_gpu(){
    if (!this -> uses_GPU()) return this;

    using namespace GPU;
    size_t N = this -> bodies.positions.size();
    clFinish(queue);
    Check(clEnqueueReadBuffer(queue, this -> positions_cl, CL_FALSE, 0, sizeof(vec3) * N, this -> bodies.positions.data(), 0, nullptr, nullptr));
    Check(clEnqueueReadBuffer(queue, this -> velocities_cl, CL_FALSE, 0, sizeof(vec3) * N, this ->bodies.velocities.data(), 0, nullptr, nullptr));

    return this;
}

Simulation::~Simulation(){
    using namespace GPU;
    clFinish(queue);
    if (!has_loaded) return;
    delete solver;
    delete updater;
    if (!uses_GPU())  return;
    for (auto k : kernels)
        clReleaseKernel(k.second);

    reset_GPU();

}

Simulation* Simulation::load(){
    if (has_loaded) return this;
    loading_queue.load(this);
    has_loaded = true;
    this -> calculate_total_mass();

    this -> GPU_on = (gravity_solver == SIMPLE_GPU || gravity_solver == GPU_BARNES_HUT);
    if (this -> GPU_on) this -> load_gpu_buffers();

    this -> updater = new Updater(this, integrator_order, integrator);
    
    this -> solver = bake_solver(this, gravity_solver);
    this -> solver -> set_force(this -> force_used);

    this -> solver -> load();

    return this;
}


Simulation* Simulation::get_force_with(force_func p){
    this -> solver -> set_force(p);
    this -> force_used = CUSTOM_F;

    return this;
}

Simulation* Simulation::update_with(Updater* p){
    if (critical_if(!p, "The given updater is null")) return this;
    this -> updater = p;
    return this;
}


Simulation* Simulation::add_graph(typename Graphs::Graph::collector_fptr listener, graph_type type){
    assert(type == BETWEEN && "cannot assign this specific function to anything other then a BETWEEN graph");
    this -> between_graphs.push_back({this, listener, type});
    this -> between_graphs.back().name = "Graph" + std::to_string(int(this -> graphs.size()));
    return this;
}

Simulation* Simulation::add_graph(typename Graphs::Graph::listener_fptr listener, graph_type type){
    assert(type != BETWEEN && "cannot assign this specific function to be BETWEEN graph");
    this -> graphs.push_back({this, listener, type});
    this -> graphs.back().name = "Graph" + std::to_string(int(this -> graphs.size()));
    return this;
}





REAL Simulation::get_height() const{
    return this -> data.size.y * get_scale();
}


bool Simulation::uses_GPU() const{
    return this -> GPU_on;
}


REAL Simulation::get_width() const{
    return this -> data.size.x * get_scale();
}


REAL Simulation::get_depth() const{
    return this -> data.size.z * get_scale();
}



REAL Simulation::get_render_height() const{
    return this -> data.size.y;
}


REAL Simulation::get_render_width() const{
    return this -> data.size.x;
}


REAL Simulation::get_render_depth() const{
    return this -> data.size.z;
}


int Simulation::get_cores() const{
    return this -> data.NUM_THREADS;
};


REAL Simulation::get_G() const{
    return this -> data.G;
};



REAL Simulation::get_c() const{
    return this -> data.c;
}


REAL Simulation::get_c_sqr() const{
    return this -> data.c_squared;
}


REAL Simulation::get_dt() const{
    return this -> data.dt;
}


REAL Simulation::get_e_sqr() const{
    return this -> data.e_squared;
}


REAL Simulation::get_takeover() const{
    return this -> data.takeover;
}


REAL Simulation::get_scale() const{
    return this -> data.simulation_scale;
}




REAL Simulation::get_heat_capacity() const{
    return this -> data.avr_heat_capacity;
}


REAL Simulation::get_boltzmann() const{
    return this -> data.boltzmann;
}


void Simulation::trigger_all_graphs(){
    for (auto& graph : graphs)
        graph.trigger();
}


bool Simulation::has_loaded_yet() const {
    return has_loaded;
}

REAL Simulation::get_total_mass() const{
    return total_mass;
}


REAL Simulation::get_adaptive_coeff() const{
    return data.adaptive_coeff;
}


Simulation* Simulation::set_adaptive_coeff(REAL s){
    data.adaptive_coeff = s;
    return this;
}



Simulation* Simulation::calculate_total_mass(){
    total_mass = 0;

    for (const auto& num : this -> bodies.masses)
        total_mass += num;
    
    return this;
}

inline vec3 Simulation::get_corner(int n) const{
    assert(n >= 0 && n < 9);

    return {
        data.size.x * (n % 2),
        data.size.y * (n==2 || n==3 || n==6 || n==7),
        data.size.z * (n > 3)
    };
}


REAL Simulation::get_time_passed() const {
    return time_passed;
}


Simulation* Simulation::collect_hamiltonian(){
    this -> add_graph(Graphs::hamiltonian_collector, ONCE);
    return this;
}


Simulation* Simulation::collect_error(){
    this -> add_graph(Graphs::error_collector, ONCE);
    return this;
}


Simulation* Simulation::collect_distance(){
    this -> add_graph(Graphs::distance_collector, FOR_EACH);
    return this;
}


Simulation* Simulation::get_force_with(force_type t){
    this -> solver -> set_force(t);
    this -> force_used = t;
    return this;
}


Simulation* Simulation::update_with(update_type t){
    this -> data.selected_update = t;
    this -> updater = new Updater(this, 0, t);
    this -> integrator = t;
    return this;
}


/*                   //===---------------------------------------------------------===//
.                    // SINGLE THREAD IMPLEMENTATION                                  //
.                    //===---------------------------------------------------------===*/


SingleThread::SingleThread(Simulation* _s){
    this -> type = SINGLE_THREAD;
    this -> _s = _s;
    this -> get_force = newtonian;
}


void SingleThread::compute_forces(){
    size_t N = this -> _s -> bodies.positions.size();
    for (int i = this->get_lower_bound(); i < this->get_upper_bound(); ++i){ 
        // saves a reference to mass, velocity and position
        double m1 = this -> _s -> bodies.masses[i];
        vec3 v1 = this -> _s -> bodies.velocities[i]; 
        vec3 p1 = this -> _s -> bodies.positions[i];

        for (int j = i+1; j < this->get_upper_bound(); ++j) {
            if (i == j) continue;

            vec3 a = this -> get_force(
                m1, this -> _s -> bodies.masses[j],
                v1, this -> _s -> bodies.velocities[j],
                p1, this -> _s -> bodies.positions[j],
                this -> _s
            );

            
            this -> _s -> bodies.accs[i] += a / m1;
            this -> _s -> bodies.accs[j] += -a / this -> _s -> bodies.masses[j];
     
        }
    }
}


/*                   //===---------------------------------------------------------===//
.                    // MULTI THREAD IMPLEMENTATION                                   //
.                    //===---------------------------------------------------------===*/


Parallelized::Parallelized(Simulation* _s){
    this -> type = PARALLEL;
    this -> _s = _s;
    this -> get_force = newtonian;
}


void Parallelized::compute_forces(){
    size_t N = this -> _s -> bodies.positions.size();
    static auto* s = this -> _s;
    auto&& get_force = this -> get_force;

    N = N - this -> lower_int_bound - ((this -> upper_int_bound < 0) ? 0 : N - this -> upper_int_bound);
    if (warn_if(N <= 0, "either upper or lower bound were set too high in the gravity solver, resulting in a negative amount of bodies to load")) 
        return; 

    N = this -> _s -> bodies.positions.size();

    parallel_for(size_t(0), size_t(N), [&](size_t i ){
        vec3 acc = vec3(0);
        double m1 = s -> bodies.masses[i];
        vec3 v1 = s -> bodies.velocities[i], p1 = s -> bodies.positions[i];
        for (int j = this -> lower_int_bound; j < ((this -> upper_int_bound < 0) ? N : this -> upper_int_bound) ; ++j) {
            if (i == j) continue;

            acc +=  this -> get_force(
                m1, s -> bodies.masses[j],
                v1, s -> bodies.velocities[j],
                p1, s -> bodies.positions[j],
                s
            ) / m1;
        }

        s -> bodies.accs[i] = acc;
    });
}
}