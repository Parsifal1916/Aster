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

template <typename T> class Barnes_Hut;

template <typename T> T pick_newpos(Simulation<T>* _s);
template <typename T> void compute_tidal_heating(Barnes_Hut<T>* _s, Body<T>* body);

template <typename T>
struct Node{
    double size = 0;
    T center_of_mass;
    T velocity;
    T center;

    size_t child = 0;
    float temp = 0;
    double mass = 0;   

    Node(T c, double s);

    void merge(Body<T>* _b);
    void add_twob(Body<T>* _b1, Node<T>* _b2);
    void init(Body<T>* _b);
    void init(Node<T>* _n);
    bool is_leaf() const;
    bool is_empty() const;
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

    void step() override;
    void update_forces() override;

    virtual void get_node_body(size_t node, Body<T>* body);

    protected:
    size_t num_childs;
    std::vector<size_t> sections;

    int opt_position(T p, T center);
    int get_to_best_leaf(Body<T>* _b);
    void insert(Body<T>* body);
    void calculate_com();

    void make_tree();
    void init_node(Node<T>& _n, Body<T>* b) const;
    void make_sections();
    size_t subdivide(int n);
    friend void compute_tidal_heating<T>(Barnes_Hut<T>* _s, Body<T>* body);
};

template <typename T>
void update_bundle(Barnes_Hut<T>* _s, unsigned short index);
}
}

#include "Aster/impl/barnes_hut.tpp"