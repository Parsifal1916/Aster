#pragma once
#include <vector>
#include <random>
#include <map>
#include <thread>

#include "Aster/simulations/basic.h"

namespace Aster{

//used to generate a random angle (0-360)
extern std::uniform_real_distribution<double> angle_rnd;

// returns a double from 0 to 1
extern std::uniform_real_distribution<double> normalized_rnd;

/*
* basic gravity solver using a single thread
*/
template <typename T>
class SingleThread final : public Simulation<T> {
    public:
    // generate from a predeterminated meta obj
    SingleThread(sim_meta m);
    SingleThread();

    void step() override;
    void update_forces() override;
};

template <typename T>
class Parallelized final : public Simulation<T> {
    public:
    // contains the necessary threads, it is reserved for _s -> get_cores() items
    std::vector<std::thread> threads;

    Parallelized(sim_meta m);
    Parallelized();
    
    void step() override;
    void update_forces() override;
    
    /*
    * sets the maximum number of threads for the given simulation
    * if the given number is higher than the number of objs
    * it will be topped at the number of objs
    */
    Parallelized<T>* set_max_threads(unsigned int t_);
};

}
