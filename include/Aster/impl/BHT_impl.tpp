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
    this -> update_bodies = get_update_func<T>(this -> data.selected_update);
    this -> data.graph_height *= this -> get_height();

    this -> threads.reserve(this -> get_cores());
    this -> obj = this -> bodies.positions.size(); 
}


//===---------------------------------------------------------===//
// Inserting body into tree                                      //
//===---------------------------------------------------------===//

template <typename T>
void BHT<T>::get_node_body(size_t node, size_t index){

    assert(node < this -> nodes.size() - 1);
    auto& cnode = this -> nodes[node];
    
    if (cnode.is_empty()) return;

    double d_squared = (this -> bodies.get_position_of(index) - this -> nodes[node].center_of_mass).sqr_magn();
    
    if (cnode.is_leaf() || d_squared * this -> theta * this -> theta > this -> nodes[node].size * this -> nodes[node].size){

        if (this -> nodes[node].center_of_mass == this -> bodies.get_position_of(index) || !cnode.mass) return;

        this -> bodies.get_acc_of(index) += this -> get_force(
            this -> bodies.get_mass_of(index),     this -> nodes[node].mass,
            this -> bodies.get_velocity_of(index), this -> nodes[node].velocity,
            this -> bodies.get_position_of(index), this -> nodes[node].center_of_mass,
            this
        ) / this -> bodies.get_mass_of(index);

        return;
    }
        
    for (int i = 0; i < this -> num_childs; i++)
        get_node_body(this -> nodes[node].child + i, index);

}
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