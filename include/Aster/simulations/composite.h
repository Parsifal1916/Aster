#pragma once
#include <vector>
#include <cassert>
#include <thread>

#include "Aster/simulations/sim_obj.h"

#include "Aster/physics/vectors.h"
#include "Aster/physics/body.h"

namespace Aster{   

template <typename T> class CompositeBody;

class Composite2d: public Simulation<vec2>{
    public:
    std::vector<std::thread> threads;
    std::vector<CompositeBody<vec2>> composites;    
    Composite2d();

    void step() override;

};

}
