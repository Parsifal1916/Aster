#pragma once
#include <vector>
#include <cassert>
#include <thread>

#include "Aster/simulations/barnes-hut.h"

namespace Aster{   
namespace Barnes{

/*                   //===---------------------------------------------------------===//
.                    // Barnes-Hut termal definition                                  //
.                    //===---------------------------------------------------------===*/



class BHT final : public Barnes_Hut {
    public: 

    BHT(sim_meta m);
    BHT();

    
};

}

}