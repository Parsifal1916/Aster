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
class BHT final : public Barnes_Hut<T> { public:
    double theta = 0.8;
    std::vector<std::thread> threads;
    std::vector<Node<T>> nodes;
    static constexpr double kelvin_to_joule = 1e-22/7.2419716667;

    BHT(sim_meta m);
    BHT();

    void step() override;

    protected:
    int num_childs;
    int opt_position(T p, T center);
    int get_to_best_leaf(Body<T>* _b);
    void insert(Body<T>* body);
    void calculate_com();
    void get_node_body(size_t node, Body<T>* body) override ;
    void update_bodies();
    void make_tree();
    void init_node(Node<T>& _n, Body<T>* b) const;
    size_t subdivide(int n);
    static double update_bundle(BHT<T>* _s, unsigned short index);
};

}
}

#include "Aster/impl/BHT_impl.tpp"