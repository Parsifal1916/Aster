#include <vector>
#include <string>
#include <fstream> 
#include <functional>

#include "Aster/graphs/graph_collection.h"
#include "Aster/simulations/sim_obj.h"

namespace Aster{

namespace Graphs{

template <typename T>
double get_total_energy(Simulation<T>* _s){
    double retval = 0.0;
    std::mutex retval_mutex;
    const int num_threads = 16;
    int num_bodies = _s->bodies.positions.size();
    
    auto thread_func = [&retval, &retval_mutex, _s, num_bodies](int start_idx, int end_idx) {
        double partial_retval = 0.0;

        for (int i = start_idx; i < end_idx; ++i) {
            partial_retval += 0.5 * _s -> bodies.get_velocity_of(i).sqr_magn() * _s -> bodies.get_mass_of(i);
            for (int j = i + 1; j < num_bodies; ++j) {
                double r = (_s -> bodies.get_position_of(j) - _s -> bodies.get_position_of(i)).magnitude();
                partial_retval -= _s->get_G() * _s -> bodies.get_mass_of(j) * _s -> bodies.get_mass_of(i) / r;
            }
        }

        std::lock_guard<std::mutex> lock(retval_mutex);
        retval += partial_retval;
    };

    std::vector<std::thread> threads;
    int chunk_size = num_bodies / num_threads;

    for (int i = 0; i < num_threads; ++i) {
        int start_idx = i * chunk_size;
        int end_idx = (i == num_threads - 1) ? num_bodies : (i + 1) * chunk_size;
        threads.push_back(std::thread(thread_func, start_idx, end_idx));
    }

    for (auto& t : threads) {
        t.join();
    }

    return retval;
}

template <typename T>
double hamiltonian_collector(Graph<T>* g, Simulation<T>* _s, size_t b){
    return get_total_energy(_s);
}

template <typename T>
double error_collector(Graph<T>* g, Simulation<T>* _s, size_t b){
    static double initial_energy = std::nan("");

    double current_energy = get_total_energy(_s);

    if (std::isnan(initial_energy)){
        initial_energy = current_energy;
        return 0;
    }

    return std::abs(initial_energy - current_energy) / std::abs(initial_energy);
}

template <typename T>
double distance_collector(Graph<T>* g, Simulation<T>* _s, size_t b){
    static std::mutex mtx;
    static std::pair<double, T> baricenter = {-1, T(0)};

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