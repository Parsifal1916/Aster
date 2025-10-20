#pragma once
#include <vector>
#include <cassert>
#include <thread>


#include "Aster/simulations/sim_obj.h"

namespace Aster{   
namespace Barnes{

/*                   //===---------------------------------------------------------===//
.                    // Nodes definition                                              //
.                    //===---------------------------------------------------------===*/

#define PRECISION_BITS 10

 class Barnes_Hut;

 vec3 pick_newpos(Simulation* _s);

struct NodesArray{
    std::vector<vec3> centers_of_mass; // center of mass of the node
    std::vector<vec3> velocities; // avrg velocity of the node
    
    std::vector<signed int> left_nodes, right_nodes;
    std::vector<REAL> temps; // avrg temperature of the node
    std::vector<REAL> masses;  // combined mass of the node

    NodesArray(){}

    /**
    * @brief merges two bodies in one node in case of conflict
    * @param body: ptr to the body to merge
    */
    void merge(Simulation* _s, size_t  _b, size_t node);

    void clear();

    void push_back(Simulation* _s, size_t index);

    /**
    * @brief allocates the right amount of memory to the nodes
    */
    void allocate(size_t n);

    size_t size() const;

    /**
    @brief resizes all arrays
    */
    void resize(size_t n);

    /**
    * @brief merges two nodes in one node
    */
    void unite(size_t start);


    /**
    * @brief initializes the node with the given body
    * @param _b: ptr to the body to copy
    */
    void init(Simulation* _s , size_t  _b, size_t index);

    /**
    * @brief returns if the node is a leaf
    * @returns child == 0
    */
    bool is_leaf(size_t i) const;

    /**
    * @brief returns if the node is empty
    * @returns mass == 0 && is_leaf()
    */
    bool is_empty(size_t i) const;
};


/*                   //===---------------------------------------------------------===//
.                    // Barnes-Hut definition                                         //
.                    //===---------------------------------------------------------===*/


class Barnes_Hut: public Solver{
    public:
    long int compressed_mortons_size  = 0;
    vec3 bounding_box;
    std::vector<std::thread> threads;
    NodesArray nodes;
    std::vector<std::pair<uint32_t, size_t>> mortons;

    Barnes_Hut(Simulation* _s);
    Barnes_Hut() {};

    Barnes_Hut* set_theta(REAL _t);

    /**
    * @brief calculates the force acting between a node and a body **recursive**+
    * @param node: node to calculate the force from
    * @param body: ptr to the body being acted upon
    * @returns nothing everything is done internally
    */
    virtual void get_node_body(signed long node, size_t  body, REAL size = 0);

    void make_tree();  
    void compute_forces() override;

    protected:

    size_t num_childs;
    std::vector<size_t> sections;



};

 
void compare_bounding_vectors(const vec3& to_compare, vec3& res);

/**
* @brief generates the morton codes for a given point
*/
 
uint32_t get_morton(Barnes_Hut* _s, vec3 point);


/**
* @brief interleaps the bits of a given number
*/
uint32_t interleap_coord(uint16_t num);

}


}