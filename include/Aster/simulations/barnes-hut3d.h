#pragma once
#include <vector>
#include <cassert>
#include <mutex>
#include <thread>

#include "Aster/physics/vectors.h"
#include "Aster/physics/body.h"
#include "Aster/physics/tool-chain.h"

#include "Aster/simulations/sim_obj.h"

#include "Aster/building-api/3d_builder.h"

namespace Aster{   
namespace Barnes{

struct Node3d{
    double size = 0;

    vec3 center_of_mass = {0,0,0};
    vec3 velocity = {0,0,0};
    vec3 center = {0,0,0};

    size_t child = 0;
    float temp = 0;
    double mass = 0;   

    Node3d(vec3 c, double s);

    void merge(Body3d* _b);
    void add_twob(Body3d* _b1, Node3d* _b2);
    void init(Body3d* _b);
    void init(Node3d* _n);
    bool is_leaf() const;
    bool is_empty() const;
    operator bool() const;
};

/*                   //===---------------------------------------------------------===//
.                    // Barnes-Hut definition                                         //
.                    //===---------------------------------------------------------===*/

class Barnes_Hut3d : public Simulation3d{public:
    double theta = 0.8;
    std::vector<std::thread> threads;
    std::vector<int> sections;
    std::vector<struct Node3d> nodes3d;

    Barnes_Hut3d(sim3d_meta m);
    Barnes_Hut3d();
    void step() override;

    private:
    int opt_position(vec3 p, vec3 center);
    int get_to_best_leaf(Body3d* _b);
    void insert(Body3d* body);
    void calculate_com();
    virtual void get_node_body(size_t node, Body3d* body);
    void update_bodies();
    void make_tree();
    void init_node(Node3d& _n, Body3d* b) const;
    size_t subdivide(int n);
    bool is_inbound(vec3& v);
    void make_sections();
    vec3 pick_newpos();
    static void update_bundle(Barnes_Hut3d* _s, unsigned short index);
};


}
}
