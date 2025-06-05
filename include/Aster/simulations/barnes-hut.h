#pragma once
#include <vector>
#include <cassert>
#include <thread>

#include "Aster/simulations/sim_obj.h"

#include "Aster/physics/vectors.h"
#include "Aster/physics/body.h"

namespace Aster{   
namespace Barnes{

/*                   //===---------------------------------------------------------===//
.                    // Nodes definition                                              //
.                    //===---------------------------------------------------------===*/

#define PRECISION_BITS 10
#define IGNORED_BITS 2

template <typename T> class Barnes_Hut;

template <typename T> T pick_newpos(Simulation<T>* _s);
template <typename T> void compute_tidal_heating(Barnes_Hut<T>* _s, size_t  body);

template <typename T>
struct NodesArray{
    std::vector<T> centers_of_mass; // center of mass of the node
    std::vector<T> velocities; // avrg velocity of the node
    
    std::vector<signed long int> left_nodes, right_nodes;
    std::vector<float> temps; // avrg temperature of the node
    std::vector<double> masses;  // combined mass of the node

    NodesArray(){}

    /**
    * @brief merges two bodies in one node in case of conflict
    * @param body: ptr to the body to merge
    */
    void merge(Simulation<T>* _s, size_t  _b, size_t node);

    void clear();

    void push_back(Simulation<T>* _s, size_t index);

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
    void init(Simulation<T>* _s , size_t  _b, size_t index);

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

template <typename T>
class Barnes_Hut: public Simulation<T>{
    public:
    double theta = 0.8;
    long int compressed_mortons_size  = 0;
    T bounding_box;
    std::vector<std::thread> threads;
    NodesArray<T> nodes;
    std::vector<std::pair<uint32_t, size_t>> mortons;

    Barnes_Hut();

    /**
    * @brief steps the simulation forward
    */
    void step() override;

    Barnes_Hut<T>* set_theta(double _t);

    /**
    * @brief calculates the force acting between a node and a body **recursive**+
    * @param node: node to calculate the force from
    * @param body: ptr to the body being acted upon
    * @returns nothing everything is done internally
    */
    virtual void get_node_body(signed long node, size_t  body, double size = 0);

    void make_tree();

    /**
    * @brief generates an array for the bodies slicing
    */
    void make_sections();
    
    Simulation<T>* use_GPU() override;

    friend void compute_tidal_heating<T>(Barnes_Hut<T>* _s, size_t  body);

    protected:

    size_t num_childs;
    std::vector<size_t> sections;



};

/**
* @brief internal function to update a batch of bodies
* @param _s: ptr to the simulation
* @param index: what chunck of the bodies array to scan
*/
template <typename T>
void update_bundle(Barnes_Hut<T>* _s, unsigned short index);

template <typename T> 
void compare_bounding_vectors(const T& to_compare, T& res);

/**
* @brief generates the morton codes for a given point
*/
template <typename T> 
uint32_t get_morton(Barnes_Hut<T>* _s, T point);


/**
* @brief interleaps the bits of a given number
*/
uint32_t interleap_coord(uint16_t num, size_t shift = 1);


template <typename T>
void barnes_update_forces(Simulation<T>* _s);


}


}

#include "Aster/impl/barnes_hut.tpp"