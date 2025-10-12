#include <vector>
#include <string>
#include <fstream> 
#include <functional>

#include "Aster/graphs/graph_collection.h"
#include "Aster/simulations/sim_obj.h"

#include <tbb/parallel_for.h>
using namespace tbb;

namespace Aster{

namespace Graphs{


REAL get_total_energy(Simulation* _s){
    REAL retval = 0.0;
    std::mutex retval_mutex;
    const int num_threads = 16;
    int num_bodies = _s->bodies.positions.size();
    
    parallel_for(size_t(0),  size_t(num_bodies), [&retval, &retval_mutex, _s, num_bodies](size_t i) {
        REAL partial_retval = 0.5 * _s -> bodies.get_velocity_of(i).sqr_magn() * _s -> bodies.get_mass_of(i);
        for (int j = i + 1; j < num_bodies; ++j) {
            REAL r = (_s -> bodies.get_position_of(j) - _s -> bodies.get_position_of(i)).magnitude();
            partial_retval -= _s->get_G() * _s -> bodies.get_mass_of(j) * _s -> bodies.get_mass_of(i) / (r+.1);
        }
        
        std::lock_guard<std::mutex> lock(retval_mutex);
        retval += partial_retval;
    });

    return retval;
}


REAL hamiltonian_collector(Graph* g, Simulation* _s, size_t b){
    return get_total_energy(_s);
}


REAL error_collector(Graph* g, Simulation* _s, size_t b){
    static REAL initial_energy = std::nan("");

    REAL current_energy = get_total_energy(_s);

    if (std::isnan(initial_energy)){
        initial_energy = current_energy;
        return 0;
    }

    return std::abs(initial_energy - current_energy) / std::abs(initial_energy);
}


REAL distance_collector(Graph* g, Simulation* _s, size_t b){
    static std::mutex mtx;
    static std::pair<REAL, vec3> baricenter = {-1, vec3(0)};

    if (baricenter.first != _s -> get_time_passed()){
        mtx.lock();
        if (baricenter.first != _s -> get_time_passed()){
            baricenter.first = _s -> get_time_passed();
            baricenter.second = vec3(0);

            for (int i = 0; i < _s -> bodies.positions.size(); ++i)
                baricenter.second += _s -> bodies.get_position_of(i) * _s -> bodies.get_mass_of(i);

            baricenter.second = baricenter.second / _s -> get_total_mass();
        }
        mtx.unlock();
    }

    return (baricenter.second - _s -> bodies.get_position_of(b) + .1).magnitude();
}

}
}