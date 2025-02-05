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
    this -> get_force = force_funcs<T>[this -> data.selected_force];
    this -> update_body = update_funcs<T>[this -> data.selected_update];
    this -> data.graph_height *= this -> get_height();

    this -> threads.reserve(this -> get_cores());
    this -> obj = this -> bodies.size(); 
}


//===---------------------------------------------------------===//
// Inserting body into tree                                      //
//===---------------------------------------------------------===//

template <typename T>
void BHT<T>::get_node_body(size_t node, Body<T>* body){

    assert(node < this -> nodes.size() - 1);
    auto& cnode = this -> nodes[node];
    
    if (cnode.is_empty()) return;

    double d_squared = (body -> position - this -> nodes[node].center_of_mass).sqr_magn();
    
    if (cnode.is_leaf() || d_squared * this -> theta * this -> theta > this -> nodes[node].size * this -> nodes[node].size){

        if (this -> nodes[node].center_of_mass == body -> position || !cnode.mass) return;

        body -> acceleration += this -> get_force(
            body -> mass, this -> nodes[node].mass,
            body -> velocity, this -> nodes[node].velocity,
            body -> position, this -> nodes[node].center_of_mass,
            this
        ) / body -> mass;

        return;
    }
        
    for (int i = 0; i < this -> num_childs; i++)
        get_node_body(this -> nodes[node].child + i, body);

}
template <typename T>
void update_bundle(BHT<T>* _s, unsigned short index){
    unsigned int mult, start, stop;
    mult = _s -> obj/_s -> get_cores();
    start = index * mult;
    stop = (index + 1) * mult;

    stop = (stop + mult > _s -> obj) ? _s -> obj : stop;
              
    for (int i = start; i < stop; ++i){  
        Body<T>* body = &_s -> bodies[i];
        body -> acceleration.reset();
        _s -> get_node_body(0, body);
        _s -> update_body(body, _s);
        
        compute_tidal_heating(_s, body);
    }
}

//

}
}