#include <vector>
#include <string>
#include <fstream> 
#include <functional>

#include "Aster/graphs/graph_collection.h"
#include "Aster/simulations/sim_obj.h"


#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
using namespace tbb;

namespace Aster{

namespace Graphs{

REAL get_total_energy(Simulation* _s) {
    const int num_bodies = _s->bodies.positions.size();
    const REAL G = _s->get_G();

    std::vector<vec3> positions(num_bodies);
    std::vector<vec3> velocities(num_bodies);
    std::vector<REAL> masses(num_bodies);

    for (int i = 0; i < num_bodies; ++i) {
        positions[i] = _s->bodies.get_position_of(i);
        velocities[i] = _s->bodies.get_velocity_of(i);
        masses[i]     = _s->bodies.get_mass_of(i);
    }

    REAL total_energy = tbb::parallel_reduce(
        tbb::blocked_range<size_t>(0, num_bodies),
        (REAL)0.0,
        [&](const tbb::blocked_range<size_t>& range, REAL local_sum) -> REAL {
            for (size_t i = range.begin(); i < range.end(); ++i) {
                local_sum += 0.5 * masses[i] * velocities[i].sqr_magn();

                for (size_t j = i + 1; j < num_bodies; ++j) {
                    vec3 dr = positions[j] - positions[i];
                    REAL r2 = dr.sqr_magn();
                    REAL inv_r = 1.0 / (sqrt(r2) + _s -> softening);
                    local_sum -= G * masses[i] * masses[j] * inv_r;
                }
            }
            return local_sum;
        },
        std::plus<REAL>()
    );

    return total_energy;
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
    vec3 baricenter = _s->get_center_of_mass();

    return (baricenter- _s -> bodies.get_position_of(b) + .1).magnitude();
}

}
}