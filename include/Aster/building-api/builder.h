#pragma once
#include <vector>
#include <random>
#include <map>
#include <thread>

#include "Aster/simulations/basic.h"

namespace Aster{

//used to generate a random angle (0-360)
extern std::uniform_real_distribution<REAL> angle_rnd;

// returns a REAL from 0 to 1
extern std::uniform_real_distribution<REAL> normalized_rnd;

/*
* basic gravity solver using a single thread
*/
class SingleThread : public Solver {
    public:
    // generate from a predeterminated meta obj
    SingleThread(Simulation* _s);
    void load() override {};
    void compute_forces() override;
};

class Parallelized : public Solver {
    public:

    Parallelized(Simulation* _s);
    void load() override {};
    void compute_forces() override;
};

class SimpleGPU : public Solver {
    public:
    SimpleGPU(Simulation* s);
    void load() override;
    void compute_forces() override;

    private:
    func_ptr load_force_kernel;
};

}
