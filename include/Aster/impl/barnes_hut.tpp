#pragma once

#include <vector>
#include <cassert>
#include <thread>
#include <type_traits>

#include "Aster/simulations/sim_obj.h"
#include "Aster/building-api/builder.h"

#include "Aster/physics/vectors.h"
#include "Aster/physics/body.h"
#include "Aster/physics/tool-chain.h"

#include "Aster/simulations/barnes-hut.h"

namespace Aster{   
namespace Barnes{


/*                   //===---------------------------------------------------------===//
.                    // NODE IMPLEMENTATION                                           //
.                    //===---------------------------------------------------------===*/


template <typename T>
Node<T>::Node(T c, double s){
    this -> center = c;
    this -> size = s;
}

template <typename T>
bool Node<T>::is_leaf() const {
    return child == 0;
}

template <typename T>
bool Node<T>::is_empty() const {
    return !mass && is_leaf();
}

template <typename T>
void Node<T>::init(Body<T>* _b){
    assert(this -> is_empty());
    mass = _b -> mass;
    center_of_mass = _b-> position;
    velocity = _b -> velocity;
    temp = _b -> temp;
}

template <typename T>
void Node<T>::init(Node<T>* _n){
    assert(this -> is_empty());
    mass = _n -> mass;
    center_of_mass = _n-> center_of_mass;
    velocity = _n -> velocity;
    temp = _n -> temp;
}

template <typename T>
void Node<T>::merge(Body<T>* _b){
    assert(child == 0);
    center_of_mass = (_b -> position * _b -> mass + center_of_mass * mass) / (mass + _b -> mass); 
    mass += _b -> mass;
    velocity += _b -> velocity;
    temp = (temp + _b -> temp) / 2;
}

template <typename T>
void Node<T>::add_twob(Body<T>* _b, Node<T>* _n){
    assert(this -> is_empty());
    mass = _b -> mass + _n -> mass;
    velocity = _b -> velocity + _n -> velocity;
    temp = (_b -> temp + _n -> temp)/2;
    center_of_mass = (_b -> position * _b -> mass  + _n -> center_of_mass * _n -> mass) / mass;
}

/*                   //===---------------------------------------------------------===//
.                    // BARNES HUT IMPLEMENTATION                                     //
.                    //===---------------------------------------------------------===*/

template <typename T>
Barnes_Hut<T>::Barnes_Hut(){
    static_assert(std::is_same<T, vec2>::value || std::is_same<T, vec3>::value, "Invalid type for class construction");

    num_childs = 4;

    if constexpr (std::is_same<T, vec3>::value)
        num_childs = 8;
    
    this -> data = sim_meta();
    this -> get_force = force_funcs<T>[this -> data.selected_force];
    this -> update_body = update_funcs<T>[this -> data.selected_update];
    this -> data.graph_height *= this -> data.size.y;

    this -> threads.reserve(this -> data.NUM_THREADS);
    this -> obj = this -> bodies.size(); 
}

template <typename T>
void Barnes_Hut<T>::step(){
    make_sections();
    make_tree();
    calculate_com();
        
    update_bodies();
    nodes.clear();

    this -> trigger_all_graphs();
    this -> time_passed++;
}

template <typename T> 
void Barnes_Hut<T>::make_sections(){
    sections.clear();
    int mult = this -> bodies.size() / this -> data.NUM_THREADS;

    for (int i = 0; i < this -> data.NUM_THREADS-1; ++i)
        sections.emplace_back(i*mult);

    sections.emplace_back(this -> bodies.size());
}

//===---------------------------------------------------------===//
// Inserting body into tree                                      //
//===---------------------------------------------------------===//

template <typename T>
void Barnes_Hut<T>::insert(Body<T>* body){
    int node_index = get_to_best_leaf(body);
    size_t cnode = node_index;

    if (nodes[cnode].is_empty()){
        nodes[cnode].init(body); //* it gets initialized
        return;
    }
    //* from now on we know it is not empty therefore it has a body ptr != nullptr

    if ((body -> position - nodes[cnode].center_of_mass).sqr_magn() < 100){ // are they in the same place?
        nodes[cnode].merge(body);
        return;
    }

    int tries = 5, new_node = node_index;
    size_t child = 0;
    short opt_1, opt_2;

    while (tries--){// trova un modo migliore
        // subdivides the current node
        child = subdivide(new_node);

        // gets where each body (the already existing one and the one we want to insert) wants to be
        opt_1 = opt_position(nodes[node_index].center_of_mass, nodes[node_index].center);
        opt_2 = opt_position(body -> position, nodes[node_index].center);

        if (opt_1 == opt_2) { // conflict
            new_node = child + opt_1;
            continue;
        }

        size_t
            p1 = child + opt_1,
            p2 = child + opt_2
        ;

        nodes[p1].init(&nodes[cnode]);
        nodes[p2].init(body);
        return;
    }

    nodes[new_node].add_twob(body, &nodes[cnode]);

}

template <typename T>
int Barnes_Hut<T>::get_to_best_leaf(Body<T>* _b){
    /*
    * NOTE: there HAS to be a leaf node at the end of the vector
    * because every time a node is created it is initialized 
    * with "is_leaf = true" and if it's a leaf it does not have childrens
    * halting the loop.
    */
    size_t node_index = 0;

    // make sure the current node is a leaf
    while (!nodes[node_index].is_leaf()) {
        short displacement = opt_position(_b -> position, nodes[node_index].center);
        node_index = nodes[node_index].child + displacement; // skips to the next optimal node
    }    

    return node_index;
}


template <typename T>
void Barnes_Hut<T>::calculate_com(){
    for (int i = this -> nodes.size()-1; i >= 0; i--){
        if (nodes[i].is_leaf()) continue;

        nodes[i].mass = 0;
        nodes[i].center_of_mass.reset();

        for (int dis = 0; dis < num_childs; ++dis){
            size_t index = nodes[i].child +dis;

            auto& child = nodes[index];

            nodes[i].mass += child.mass;
            nodes[i].center_of_mass += child.center_of_mass*child.mass;
        }

        nodes[i].center_of_mass /= nodes[i].mass;
    }
}

template <typename T>
void Barnes_Hut<T>::get_node_body(size_t node, Body<T>* body){
    assert(node < nodes.size() - 1);
    auto& cnode = nodes[node];
    
    if (cnode.is_empty()) return;

    double d_squared = (body -> position - nodes[node].center_of_mass).sqr_magn();

    if (d_squared * theta * theta > nodes[node].size * nodes[node].size){ // use the optmisation
        body -> acceleration += this -> get_force(
            body -> mass, nodes[node].mass,
            body -> velocity, nodes[node].velocity,
            body -> position, nodes[node].center_of_mass,
            this
        ) / body -> mass;

        return;
    } 
    
    if (cnode.is_leaf()){

        if (nodes[node].center_of_mass == body -> position || !cnode.mass) return;

        body -> acceleration += this -> get_force(
            body -> mass, nodes[node].mass,
            body -> velocity, nodes[node].velocity,
            body -> position, nodes[node].center_of_mass,
            this
        ) / body -> mass;
            
        return;
    }
        
    for (int i = 0; i < num_childs; i++)
        get_node_body(nodes[node].child + i, body);
    
}

template <typename T>
void update_bundle(Barnes_Hut<T>* _s, unsigned short index){
    unsigned int mult, start, stop;
    mult = _s -> obj/_s -> data.NUM_THREADS;
    start = index * mult;
    stop = (index + 1) * mult;

    stop = (stop + mult > _s -> obj) ? _s -> obj : stop;
              
    for (int i = start; i < stop; ++i){  
        Body<T>* body = &_s -> bodies[i];
        body -> acceleration.reset();
        _s -> get_node_body(0, body);
        _s -> update_body(body, _s);

        body -> temp -= body -> temp * body -> temp * body -> temp * body -> temp * _s -> data.boltzmann * _s -> data.dt;
    }
}

template <typename T>
void Barnes_Hut<T>::update_bodies(){
    this -> obj = this ->bodies.size();

    for (int i = 0; i < this -> data.NUM_THREADS; ++i)
        this -> threads.emplace_back(std::thread(update_bundle<T> , this, i));

    for (auto& t : threads)
        t.join();

    this -> threads.clear();
}


template <typename T>
void Barnes_Hut<T>::make_tree(){
    this -> threads.clear();

    nodes.push_back(Node<T>(pick_newpos<T>(this), this -> get_center().x * 2));
    
    for (auto& b : this -> bodies)
        insert(&b);
}


}
}
