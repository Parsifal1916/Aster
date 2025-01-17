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

class BHT : public Simulation {public:
    double theta = 0.8;
    std::vector<std::thread> threads;
    std::vector<struct Node> nodes;
    static constexpr double kelvin_to_joule = 1e-22/7.2419716667;

    BHT(sim_meta m);
    BHT();

    void step() override;

    private:
    int opt_position(vec2 p, vec2 center);
    int get_to_best_leaf(Body* _b);
    void insert(Body* body);
    void calculate_com();
    virtual void get_node_body(size_t node, Body* body);
    void update_bodies();
    void make_tree();
    void init_node(Node& _n, Body* b) const;
    size_t subdivide(int n);
    static double update_bundle(BHT* _s, unsigned short index);
};

}
}