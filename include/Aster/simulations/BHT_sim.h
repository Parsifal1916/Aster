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

    BHT(sim_meta m){
        this -> data = m;
        get_force = force_funcs[m.selected_force];
        update_body = update_funcs[m.selected_update];
        m.graph_height *= m.HEIGHT;

        get_rndX = std::uniform_real_distribution<double>(0, data.WIDTH);
        get_rndY = std::uniform_real_distribution<double>(0, data.HEIGHT);

        threads.reserve(data.NUM_THREADS);
        obj = bodies.size(); 
    }


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

/*                   //===---------------------------------------------------------===//
.                    // BARNES HUT IMPLEMENTATION                                     //
.                    //===---------------------------------------------------------===*/

void BHT::step(){
    make_tree();
    calculate_com();
        
    update_bodies();
    nodes.clear();
}

//===---------------------------------------------------------===//
// Inserting body into tree                                      //
//===---------------------------------------------------------===//

void BHT::insert(Body* body){
    int node_index = get_to_best_leaf(body);
    Node& cnode = nodes[node_index];

    if (cnode.is_empty()){
        cnode.init(body); //* it gets initialized
        return;
    }
    //* from now on we know it is not empty therefore it has a body ptr != nullptr

    if ((body -> position - cnode.center_of_mass).sqr_magn() < 100){ // are they in the same place?
        cnode.merge(body);
        return;
    }

    int tries = 5, new_node = node_index;
    size_t child = 0;
    short opt_1, opt_2;

    while (tries--){// trova un modo migliore
        // subdivides the current node
        child = subdivide(node_index);

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

        nodes[p1].init(&cnode);
        nodes[p2].init(body);
        return;
    }

    nodes[child].add_twob(body, &nodes[node_index]);

}

int BHT::get_to_best_leaf(Body* _b){
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


/*
* returns the index of the child node
* given the position of the object to 
* insert and the center of the root node
*
* @param p: position of the object
* @param center: center of the node
* @retval: displacement in the node vector to find the optimal node
*/
int BHT::opt_position(vec2 p, vec2 center){ 
    bool 
        x = p.x > center.x,
        y = p.y > center.y
    ;

    return x + y*2;
}

   
void BHT::calculate_com(){
    for (int i = nodes.size()-1; i >= 0; i--){
        if (nodes[i].is_leaf()) continue;

        nodes[i].mass = 0;
        nodes[i].center_of_mass = {0,0};

        for (int dis = 0; dis < 4; ++dis){
            size_t index = nodes[i].child +dis;

            auto& child = nodes[index];

            nodes[i].mass += child.mass;
            nodes[i].center_of_mass += child.center_of_mass*child.mass;
        }

        nodes[i].center_of_mass /= nodes[i].mass;
    }
}

void BHT::get_node_body(size_t node, Body* body){
    assert(node < nodes.size());
    auto& cnode = nodes[node];
    
    if (cnode.is_empty()) return;

    double
       dx = cnode.center_of_mass.x - body -> position.x,
       dy = cnode.center_of_mass.y - body -> position.y,
       d_squared = dx*dx+dy*dy
    ;

    if (d_squared * theta * theta > nodes[node].size * nodes[node].size){ // use the optmisation

        body -> acceleration += get_force(
            body -> mass, nodes[node].mass,
            body -> velocity, nodes[node].velocity,
            body -> position, nodes[node].center_of_mass,
            this
        ) / body -> mass;

        return;
    } 
    
    if (cnode.is_leaf()){

        if (nodes[node].center_of_mass == body -> position || !cnode.mass) return;

        body -> acceleration += get_force(
            body -> mass, nodes[node].mass,
            body -> velocity, nodes[node].velocity,
            body -> position, nodes[node].center_of_mass,
            this
        ) / body -> mass;


        return;
    }

    for (int i = 0; i < 4; i++)
        get_node_body(nodes[node].child + i, body);

}

void BHT::update_bodies(){
    std::vector<double> temps = {};
    temps.reserve(16);
    this -> obj = bodies.size(); 

    for (int i = 0; i < this -> data.NUM_THREADS; ++i)
        this -> threads.emplace_back([i, &temps, this] {temps.push_back(update_bundle(this, i));});

    for (auto& t : threads)
        t.join();

    this -> threads.clear();

    // i know i can use std::reduce
    // window's compiler is complaining because versions

    max_temp = 0;
    
    for (const auto& t : temps)
        max_temp = (t > max_temp) ?  t : max_temp;

    //std::cout << max_temp;
}

double BHT::update_bundle(BHT* _s, unsigned short index){
    double t = -1;
    unsigned int mult, start, stop;
    mult = _s -> obj/_s -> data.NUM_THREADS;
    start = index * mult;
    stop = (index + 1) * mult;
              
    for (int i = start; i < stop; ++i){  
        Body* body = &_s -> bodies[i];
        _s -> get_node_body(0, body);
        body -> temp = 1;//sqr_magn(body -> velocity) + 1;
        _s -> update_body(body, _s);

        
        //std::cout << magnitude(body -> acceleration) << "\n";
       //kelvin_to_joule*body -> temp * body -> temp * body -> temp * body -> temp * _s -> data.boltzmann + sqr_magn(body -> acceleration);

        if (body -> temp > t) t = body -> temp;
    }

    return t;
}

void BHT::make_tree(){
    //nodes.reserve(obj * 3);
    nodes.push_back(Node({data.WIDTH/2, data.HEIGHT/2}, data.WIDTH));
    for (auto& b : bodies){
        insert(&b);
    }
}

size_t BHT::subdivide(int n){
    //std::cout << nodes.size() << "\n";
    size_t child_index = nodes.size();
    Node& cnode = nodes[n]; 

    cnode.child = child_index;

     for (int i = 0; i < 4; i++){
        nodes.push_back(Node(
            {
            /*x = */ (0.5 - (i % 2 == 0)) * cnode.size + cnode.center.x,
            /*y = */ (0.5 - (i < 2)) * cnode.size + cnode.center.y
            },
            cnode.size/2
        ));
    }

    return child_index;
}

}
}