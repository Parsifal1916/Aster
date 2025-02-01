#pragma once
#include <vector>
#include <random>
#include <map>
#include <thread>

#include "Aster/simulations/basic.h"
#include "Aster/physics/tool-chain.h"

namespace Aster{


extern std::uniform_real_distribution<double> angle_rnd;
extern std::uniform_real_distribution<double> normalized_rnd;

template <typename T>
extern std::map<std::string, func_ptr<T>> update_funcs;

template <typename T>
extern std::map<std::string, force_func<T>> force_funcs;

template <typename T>
class SingleThread final : public Simulation<T> {
    public:
    SingleThread(sim_meta m);
    SingleThread();

    void step() override ;
};

template <typename T>
class Parallelized final : public Simulation<T> {
    public:
    std::vector<std::thread> threads;

    Parallelized(sim_meta m);
    Parallelized();
    void step() override;
    
    /*
    * sets the maximum number of threads for the given simulation
    * if the given number is higher than the number of objs
    * it will be topped at the number of objs
    */
    Parallelized<T>* set_max_threads(unsigned int t_);
};

}
