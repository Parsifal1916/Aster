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
    int num_bodies = _s->bodies.size();
    
    auto thread_func = [&retval, &retval_mutex, _s, num_bodies](int start_idx, int end_idx) {
        double partial_retval = 0.0;

        for (int i = start_idx; i < end_idx; ++i) {
            const auto& b1 = _s->bodies[i];
            partial_retval += 0.5 * b1.velocity.sqr_magn() * b1.mass;
            for (int j = i + 1; j < num_bodies; ++j) {
                const auto& b2 = _s->bodies[j];
                double r = (b2.position - b1.position).magnitude();
                partial_retval -= _s->get_G() * b1.mass * b2.mass / r;
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
double hamiltonian_collector(Graph<T>* g, Simulation<T>* _s, Body<T>* b){
    return get_total_energy(_s);
}

template <typename T>
double error_collector(Graph<T>* g, Simulation<T>* _s, Body<T>* b){
    static double initial_energy = std::nan("");

    double current_energy = get_total_energy(_s);

    if (std::isnan(initial_energy)){
        initial_energy = current_energy;
        return 0;
    }

    return std::abs(initial_energy - current_energy) / std::abs(initial_energy);
}

}
}