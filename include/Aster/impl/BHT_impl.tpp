#pragma once

#include "Aster/simulations/BHT_sim.h"

namespace Aster{
namespace Barnes{


/*                   //===---------------------------------------------------------===//
.                    // BARNES HUT TERMAL IMPLEMENTATION                              //
.                    //===---------------------------------------------------------===*/

template <typename T>
BHT<T>::BHT(){
    static_assert(std::is_same<T, vec2>::value || std::is_same<T, vec3>::value, "Invalid type for class construction");

    this -> num_childs = 4;

    if constexpr (std::is_same<T, vec3>::value)
        this -> num_childs = 8;
    
    this -> data = sim_meta();
    this -> data.type = BH_termal;
    this -> get_force = get_force_func<T>(this -> data.selected_force);
    this -> update_bodies = get_update_func<T>(this -> data.selected_update, this -> uses_GPU());
    this -> data.graph_height *= this -> get_height();

    this -> threads.reserve(this -> get_cores());
    this -> obj = this -> bodies.positions.size(); 
}


//===---------------------------------------------------------===//
// Inserting body into tree                                      //
//===---------------------------------------------------------===//

template <typename T>
void update_bundle(BHT<T>* _s, unsigned short index){
    unsigned int mult, start, stop;
    mult = _s -> obj/_s -> get_cores();
    start = index * mult;
    stop = (index + 1) * mult;

    stop = (stop + mult > _s -> obj) ? _s -> obj : stop;
              
    for (int i = start; i < stop; ++i){  
        _s -> bodies.get_acc_of(i).reset();
        _s -> get_node_body(0, i);
        _s -> update_body(i, _s);
        
        compute_tidal_heating(_s, i);
    }
}

//

}
}