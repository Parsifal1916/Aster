#pragma once

#include <vector>
#include <cassert>
#include <thread>
#include <type_traits>

#include "Aster/simulations/sim_obj.h"
#include "Aster/building-api/builder.h"

#include "Aster/physics/vectors.h"
#include "Aster/physics/body.h"

#include "Aster/simulations/barnes-hut.h"

namespace Aster{   

template <typename T> void get_new_temp(Simulation<T>* _s, size_t body, T pos, T vel, double temp, double mass);

namespace Barnes{


/*                   //===---------------------------------------------------------===//
.                    // NODE IMPLEMENTATION                                           //
.                    //===---------------------------------------------------------===*/


template <typename T>
Node<T>::Node(T c, double s){
    this -> center = c;
    this -> size = s;
}

/**
* @brief returns if the node is a leaf
* @returns child == 0
*/
template <typename T>
bool Node<T>::is_leaf() const {
    return child == 0;
}

/**
* @brief returns if the node is empty
* @returns mass == 0 && is_leaf()
*/
template <typename T>
bool Node<T>::is_empty() const {
    return !mass && is_leaf();
}

/**
* @brief initializes the node with the given body
* @param _b: ptr to the body to copy
*/
template <typename T>
void Node<T>::init(Simulation<T>* _s,size_t  _b){
    assert(this -> is_empty());
    mass           = _s -> bodies.get_mass_of(_b);
    center_of_mass = _s -> bodies.get_position_of(_b);
    velocity       = _s -> bodies.get_velocity_of(_b);
    temp           = _s -> bodies.get_temp_of(_b);
}

/**
* @brief initializes the node with the given node
* @param _n: ptr to the node to copy
*/
template <typename T>
void Node<T>::init(Node<T>* _n){
    assert(this -> is_empty());
    mass = _n -> mass;
    center_of_mass = _n-> center_of_mass;
    velocity = _n -> velocity;
    temp = _n -> temp;
}

/**
* @brief merges two bodies in one node in case of conflict
* @param body: index to the body to merge
*/
template <typename T>
void Node<T>::merge(Simulation<T>* _s, size_t _b){
    if (warn_if(child != 0, "assertion child != 0 failed in Node<T>::merge"))
        return;
    
    // combines the key values
    center_of_mass = (_s -> bodies.get_position_of(_b) * _s -> bodies.get_mass_of(_b) + center_of_mass * mass) / (mass + _s -> bodies.get_mass_of(_b)); 
    mass +=  _s -> bodies.get_mass_of(_b);
    velocity += _s -> bodies.get_velocity_of(_b);

    // for the temperature it uses the maximum
    temp = std::max(double(temp), _s -> bodies.get_temp_of(_b));
}

/**
* @brief puts a body and a node in one node
* @param _b: index to the body
* @param _n: ptr to the node
*/
template <typename T>
void Node<T>::add_twob(Simulation<T>* _s, size_t _b, Node<T>* _n){
    if (warn_if(!this -> is_empty(), "assertion child != 0 failed in Node<T>::merge"))
        return;

    // combines the values
    mass = _s -> bodies.get_mass_of(_b) + _n -> mass;
    velocity = _s -> bodies.get_velocity_of(_b) + _n -> velocity;
    temp = std::max(double(temp), _s -> bodies.get_temp_of(_b));

    // gets the combined center of mass
    center_of_mass = (_s -> bodies.get_position_of(_b)* _s -> bodies.get_mass_of(_b)  + _n -> center_of_mass * _n -> mass) / mass;
}

/*                   //===---------------------------------------------------------===//
.                    // BARNES HUT IMPLEMENTATION                                     //
.                    //===---------------------------------------------------------===*/

template <typename T>
Barnes_Hut<T>::Barnes_Hut(){
    static_assert(std::is_same<T, vec2>::value || std::is_same<T, vec3>::value, "Invalid type for class construction");

    num_childs = 4; // sets up the child number

    if constexpr (std::is_same<T, vec3>::value) // if its 3d it sets it up for an oct-tree
        num_childs = 8;
    
    // sets up the simulation
    this -> data = sim_meta(); 
    this -> data.type = BARNES_HUT;
    this -> get_force = get_force_func<T>(this -> data.selected_force);
    this -> update_bodies = get_update_func<T>(this -> data.selected_update);
    this -> data.graph_height *= this -> data.size.y;
    this -> update_forces = parallel_fu;

    this -> threads.reserve(this -> get_cores());
    this -> obj = this -> bodies.positions.size(); 
}

/**
* @brief steps the simulation forward
*/
template <typename T>
void Barnes_Hut<T>::step(){
    make_sections(); // for threading
    make_tree(); // generates the quad-tree / oct-tree
    calculate_com(); // calculates the center of mass
        
    // calculates the forces
    this -> update_forces(this);
    nodes.clear(); // cleans tree and updates bodies
    this -> update_bodies(this);

    // triggers graph and steps the time
    this -> trigger_all_graphs();
    this -> time_passed++;
}

/**
* @brief generates an array for the bodies slicing
*/
template <typename T> 
void Barnes_Hut<T>::make_sections(){
    sections.clear();
    int mult = this -> bodies.positions.size() / this -> get_cores();

    for (int i = 0; i < this -> get_cores()-1; ++i)
        sections.emplace_back(i*mult);

    sections.emplace_back(this -> bodies.positions.size());
}

//===---------------------------------------------------------===//
// Inserting body into tree                                      //
//===---------------------------------------------------------===//

/**
* @brief inserts a body into the tree
* @param body: ptr to the body to insert
*/
template <typename T>
void Barnes_Hut<T>::insert(size_t index){
    int node_index = get_to_best_leaf(index); // gets the right successor
    size_t cnode = node_index; // saves the node

    // checks if the node is empty
    if (nodes[cnode].is_empty()){
        nodes[cnode].init(this, index); // it gets initialized
        return;
    }
    // from now on we know it is not empty therefore it has a body ptr != nullptr

    if ((this -> bodies.get_position_of(index) - nodes[cnode].center_of_mass).sqr_magn() < 100){ // are they in the same place?
        nodes[cnode].merge(this, index);
        return;
    }

    int tries = 5, new_node = node_index;
    size_t child = 0;
    short opt_1, opt_2;

    while (tries--){
        // subdivides the current node
        child = subdivide(new_node);

        // gets where each body (the already existing one and the one we want to insert) wants to be
        opt_1 = opt_position(nodes[node_index].center_of_mass, nodes[node_index].center);
        opt_2 = opt_position(this -> bodies.get_position_of(index), nodes[node_index].center);

        if (opt_1 == opt_2) { // in case of conflict
            new_node = child + opt_1; // steps ahead
            continue;
        }
        
        // if the above check is false we can merge

        size_t
            p1 = child + opt_1,
            p2 = child + opt_2
        ;

        // initializes the found nodes with the respective bodies
        nodes[p1].init(&nodes[cnode]);
        nodes[p2].init(this, index);
        return;
    }

    // if the loop fails (tries <= 0) it merges them into the same node
    nodes[new_node].add_twob(this, index, &nodes[cnode]);

}

/**
* @brief finds the best index to insert the body in
* @param _b: index to the body to insert
*/
template <typename T>
int Barnes_Hut<T>::get_to_best_leaf(size_t _b){
    /*
    * NOTE: there HAS to be a leaf node at the end of the vector
    * because every time a node is created it is initialized 
    * with "is_leaf = true" and if it's a leaf it does not have childrens
    * halting the loop.
    */
    size_t node_index = 0;

    // make sure the current node is a leaf
    while (!nodes[node_index].is_leaf()) {
        short displacement = opt_position(this -> bodies.get_position_of(_b), nodes[node_index].center);
        node_index = nodes[node_index].child + displacement; // skips to the next optimal node
    }    

    return node_index;
}

/**
* @brief calculates the center of mass for all nodes
*/
template <typename T>
void Barnes_Hut<T>::calculate_com(){
    // iterates nodes backwards
    for (int i = this -> nodes.size()-1; i >= 0; i--){
        if (nodes[i].is_leaf()) continue; // skips leaf leaf

        // begins generating the COM
        nodes[i].mass = 0;
        nodes[i].center_of_mass.reset();

        for (int dis = 0; dis < num_childs; ++dis){ // adds the mass and COM of its childs 
            size_t index = nodes[i].child +dis; // updates the index
            auto& child = nodes[index];

            // calculates center of mass and adds up the masses
            nodes[i].mass += child.mass; 
            nodes[i].center_of_mass += child.center_of_mass*child.mass;
        }

        // final calculation step
        nodes[i].center_of_mass /= nodes[i].mass;
    }
}

/**
* @brief calculates the force acting between a node and a body **recursive**+
* @param node: node to calculate the force from
* @param body: ptr to the body being acted upon
* @returns nothing everything is done internally
*/
template <typename T>
void Barnes_Hut<T>::get_node_body(size_t node, size_t index){
    if (critical_if(node >= nodes.size(), "got unexpected node index in get_node_body(size_t, Body<T>*)"))
        exit(-1);

    auto& cnode = nodes[node];
    
    if (cnode.is_empty()) return;

    double d_squared = (this -> bodies.get_position_of(index) - nodes[node].center_of_mass).sqr_magn();

    if (d_squared * theta * theta > nodes[node].size * nodes[node].size){ // use the optmisation
        this -> bodies.get_acc_of(index) += this -> get_force(
            this -> bodies.get_mass_of(index),     nodes[node].mass,
            this -> bodies.get_velocity_of(index), nodes[node].velocity,
            this -> bodies.get_position_of(index), nodes[node].center_of_mass,
            this
        ) / this -> bodies.get_mass_of(index);

        return;
    } 
    
    if (cnode.is_leaf()){

        if (nodes[node].center_of_mass == this -> bodies.get_position_of(index) || !cnode.mass) return;

        this -> bodies.get_acc_of(index) += this -> get_force(
            this -> bodies.get_mass_of(index), nodes[node].mass,
            this -> bodies.get_velocity_of(index), nodes[node].velocity,
            this -> bodies.get_position_of(index), nodes[node].center_of_mass,
            this
        ) / this -> bodies.get_mass_of(index);
            
        return;
    }
        
    for (int i = 0; i < num_childs; i++)
        get_node_body(nodes[node].child + i, index);
    
}

/**
* @brief internal function to update a batch of bodies
* @param _s: ptr to the simulation
* @param index: what chunck of the bodies array to scan
*/
template <typename T>
void update_bundle(Barnes_Hut<T>* _s, unsigned short index){
    // finds the start and stop values based on the index
    unsigned int mult, start, stop;
    mult = _s -> obj/_s -> get_cores();
    start = index * mult;
    stop = (index + 1) * mult;

    // checks for an off-by-one error
    stop = (stop + mult > _s -> obj) ? _s -> obj : stop;
    
    // calculates the forces of those bodies
    for (int i = start; i < stop; ++i){  
        _s -> bodies.get_acc_of(i).reset();
        _s -> get_node_body(0, i);
    }
}

template <typename T> 
void compute_tidal_heating(Barnes_Hut<T>* _s, size_t index){
    for (int i = _s -> nodes[0].child; i < _s -> num_childs + _s -> nodes[0].child; ++i){
        if (i > _s -> bodies.positions.size()) return;
        Node<T>& cnode = _s -> nodes[i];
        
        get_new_temp<T>(_s,
            index, cnode.center_of_mass, 
            cnode.velocity, cnode.temp, 
            cnode.mass
        );
    }
}

template <typename T>
void Barnes_Hut<T>::make_tree(){
    this -> threads.clear();

    nodes.push_back(Node<T>(pick_newpos<T>(this), this -> get_center().x * 2));
    
    for (int i = 0; i < this -> bodies.positions.size(); ++i)
        insert(i);
}


}
}
