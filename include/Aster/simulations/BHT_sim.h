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


template <typename T>
class BHT final : public Barnes_Hut<T> {
    public: 

    BHT(sim_meta m);
    BHT();

    
};

template <typename T>
REAL update_bundle(BHT<T>* _s, unsigned short index);

}

}

#include "Aster/impl/BHT_impl.tpp"