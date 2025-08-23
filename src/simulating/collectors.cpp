#include <vector>
#include <string>
#include <fstream> 
#include <functional>

#include "Aster/graphs/graph_collection.h"
#include "Aster/simulations/sim_obj.h"


template <typename F>
FORCE_INLINE void parallel(int cores, size_t num, F func) {
    const size_t n_threads = std::min(num, static_cast<size_t>(cores));
    
    if (num == 0) return;
    
    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    
    const size_t chunk_size = (num + n_threads - 1) / n_threads;
    
    for (size_t i = 0; i < n_threads; ++i) {
        const size_t start = i * chunk_size;
        const size_t end = std::min(start + chunk_size, num);
        
        if (start >= end) break;
        
        threads.emplace_back([start, end, &func]() {
            for (size_t b = start; b < end; ++b) {
                func(b);
            }
        });
    }
    
    for (auto& t : threads) {
        t.join();
    }
}

namespace Aster{

namespace Graphs{

template <typename T>
REAL get_total_energy(Simulation<T>* _s){
    REAL retval = 0.0;
    std::mutex retval_mutex;
    const int num_threads = 16;
    int num_bodies = _s->bodies.positions.size();
    
    parallel(_s -> get_cores(), num_bodies, [&retval, &retval_mutex, _s, num_bodies](size_t i) {
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

template <typename T>
REAL hamiltonian_collector(Graph<T>* g, Simulation<T>* _s, size_t b){
    return get_total_energy(_s);
}

template <typename T>
REAL error_collector(Graph<T>* g, Simulation<T>* _s, size_t b){
    static REAL initial_energy = std::nan("");

    REAL current_energy = get_total_energy(_s);

    if (std::isnan(initial_energy)){
        initial_energy = current_energy;
        return 0;
    }

    return std::abs(initial_energy - current_energy) / std::abs(initial_energy);
}

template <typename T>
REAL distance_collector(Graph<T>* g, Simulation<T>* _s, size_t b){
    static std::mutex mtx;
    static std::pair<REAL, T> baricenter = {-1, T(0)};

    if (baricenter.first != _s -> get_time_passed()){
        mtx.lock();
        if (baricenter.first != _s -> get_time_passed()){
            baricenter.first = _s -> get_time_passed();
            baricenter.second = T(0);

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