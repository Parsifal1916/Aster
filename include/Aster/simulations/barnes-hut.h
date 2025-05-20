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
struct Node{
    double size = 0;
    T center_of_mass; // center of mass of the node
    T velocity; // avrg velocity of the node
    T center; // position of the node

    size_t child = 0; // index to the child in the array
    float temp = 0; // avrg temperature of the node
    double mass = 0;  // combined mass of the node

    Node(T c, double s);

    /**
    * @brief merges two bodies in one node in case of conflict
    * @param body: ptr to the body to merge
    */
    void merge(Simulation<T>* _s, size_t  _b);

    /**
    * @brief puts a body and a node in one node
    * @param _b: ptr to the body
    * @param _n: ptr to the node
    */
    void add_twob(Simulation<T>* _s, size_t _b1, Node<T>* _b2);

    /**
    * @brief initializes the node with the given body
    * @param _b: ptr to the body to copy
    */
    void init(Simulation<T>* _s , size_t  _b);

    /**
    * @brief initializes the node with the given node
    * @param _n: ptr to the node to copy
    */
    void init(Node<T>* _n);

    /**
    * @brief returns if the node is a leaf
    * @returns child == 0
    */
    bool is_leaf() const;

    /**
    * @brief returns if the node is empty
    * @returns mass == 0 && is_leaf()
    */
    bool is_empty() const;

    /**
    * @brief same as is_emtpy
    * @returns if the node is empty
    */
    operator bool() const;
};


/*                   //===---------------------------------------------------------===//
.                    // Barnes-Hut definition                                         //
.                    //===---------------------------------------------------------===*/

template <typename T>
class Barnes_Hut: public Simulation<T>{
    public:
    double theta = 0.8;
    std::vector<std::thread> threads;
    std::vector<Node<T>> nodes;

    Barnes_Hut();

    /**
    * @brief steps the simulation forward
    */
    void step() override;

    /**
    * @brief calculates the force acting between a node and a body **recursive**+
    * @param node: node to calculate the force from
    * @param body: ptr to the body being acted upon
    * @returns nothing everything is done internally
    */
    virtual void get_node_body(size_t node, size_t  body);

    protected:
    size_t num_childs;
    std::vector<size_t> sections;

    int opt_position(T p, T center);

    /**
    * @brief finds the best index to insert the body in
    * @param _b: ptr to the body to insert
    */
    int get_to_best_leaf(size_t  _b);

    /**
    * @brief inserts a body into the tree
    * @param body: ptr to the body to insert
    */
    void insert(size_t  body);

    /**
    * @brief calculates the center of mass for all nodes
    */ 
    void calculate_com();

    void make_tree();
    void init_node(Node<T>& _n, size_t  b) const;

    /**
    * @brief generates an array for the bodies slicing
    */
    void make_sections();
    size_t subdivide(int n);
    friend void compute_tidal_heating<T>(Barnes_Hut<T>* _s, size_t  body);




};

/**
* @brief internal function to update a batch of bodies
* @param _s: ptr to the simulation
* @param index: what chunck of the bodies array to scan
*/
template <typename T>
void update_bundle(Barnes_Hut<T>* _s, unsigned short index);



/**
* @brief generates the morton codes for a given point
*/
template <typename T> 
uint32_t get_morton(T point);


/**
* @brief interleaps the bits of a given number
*/
uint32_t interleap_coord(double num, size_t shift = 1);

}


}

#include "Aster/impl/barnes_hut.tpp"