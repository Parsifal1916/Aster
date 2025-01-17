#pragma once

#include <map>
#include <thread>

#include "Aster/simulations/basic.h"
#include "Aster/simulations/3d_sim_obj.h"

#include "Aster/physics/tool-chain.h"

namespace Aster{

extern std::map<std::string, func_ptr3d> update_funcs_3d;
extern std::map<std::string, force_func3d> force_funcs_3d;


class SingleThread3d : public Simulation3d {
    public:
    SingleThread3d(sim3d_meta m);
    SingleThread3d();

    void step() override;
};

class Parallelized3d : public Simulation3d {
    public:
    std::vector<std::thread> threads;
    Parallelized3d(sim3d_meta m);
    Parallelized3d();
    void step() override;
    
    /*
    * sets the maximum number of threads for the given simulation
    * if the given number is higher than the number of objs
    * it will be topped at the number of objs
    */
    Parallelized3d* set_max_threads(unsigned int t_);
    static void update_bundle(Simulation3d* _s, unsigned short index);
};

}
