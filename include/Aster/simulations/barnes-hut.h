#pragma once
#include <vector>
#include <cassert>
#include <thread>

#include "Aster/simulations/sim_obj.h"
#include "Aster/building-api/builder.h"

#include "Aster/physics/vectors.h"
#include "Aster/physics/body.h"
#include "Aster/physics/tool-chain.h"

namespace Aster{   
namespace Barnes{

/*                   //===---------------------------------------------------------===//
.                    // Nodes definition                                              //
.                    //===---------------------------------------------------------===*/

struct Node{
    double size = 0;
    vec2 center_of_mass = {0,0};
    vec2 velocity = {0,0};
    vec2 center = {0,0};

    size_t child = 0;
    float temp = 0;
    double mass = 0;   

    Node(vec2 c, double s);

    void merge(Body* _b);
    void add_twob(Body* _b1, Node* _b2);
    void init(Body* _b);
    void init(Node* _n);
    bool is_leaf() const;
    bool is_empty() const;
    operator bool() const;
};


/*                   //===---------------------------------------------------------===//
.                    // Barnes-Hut definition                                         //
.                    //===---------------------------------------------------------===*/

class Barnes_Hut : public Simulation{
    public:
    double theta = 0.8;
    std::vector<std::thread> threads;
    std::vector<struct Node> nodes;

    Barnes_Hut(sim_meta m);
    Barnes_Hut();

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
    friend void update_bundle(Barnes_Hut* _s, unsigned short index);
};

}
}
