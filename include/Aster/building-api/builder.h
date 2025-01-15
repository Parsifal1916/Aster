#pragma once
#include <vector>
#include <random>
#include <map>

#include "Aster/simulations/basic.h"
#include "Aster/physics/tool-chain.h"

namespace Aster{


extern std::uniform_real_distribution<double> angle_rnd;
extern std::uniform_real_distribution<double> normalized_rnd;

extern std::map<std::string, func_ptr> update_funcs;
extern std::map<std::string, force_func> force_funcs;


class SingleThread : public Simulation {
    public:
    SingleThread(sim_meta m);
    SingleThread();

    void step() override ;
};

class Parallelized : public Simulation {
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
    Parallelized* set_max_threads(unsigned int t_);
    static void update_bundle(Simulation* _s, unsigned short index);
};

}
