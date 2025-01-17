#include <vector>
#include <cassert>
#include <mutex>
#include <thread>

#include "Aster/physics/vectors.h"
#include "Aster/physics/body.h"
#include "Aster/physics/tool-chain.h"

#include "Aster/simulations/sim_obj.h"
#include "Aster/simulations/barnes-hut3d.h"

#include "Aster/building-api/3d_builder.h"

namespace Aster{   
namespace Barnes{

/*                   //===---------------------------------------------------------===//
.                    // NODE IMPLEMENTATION                                           //
.                    //===---------------------------------------------------------===*/

Node3d::Node3d(vec3 c, double s){
    center = c;
    size = s;
}

bool Node3d::is_leaf() const {
    return child == 0;
}

bool Node3d::is_empty() const {
    return !mass && is_leaf();
}

void Node3d::init(Body3d* _b){
    assert(this -> is_empty());
    mass = _b -> mass;
    center_of_mass = _b-> position;
    velocity = _b -> velocity;
    temp = _b -> temp;
}

void Node3d::init(Node3d* _n){
    assert(this -> is_empty());
    mass = _n -> mass;
    center_of_mass = _n-> center_of_mass;
    velocity = _n -> velocity;
    temp = _n -> temp;
}

void Node3d::merge(Body3d* _b){
    assert(child == 0);
    center_of_mass = (_b -> position * _b -> mass + center_of_mass * mass) / (mass + _b -> mass); 
    mass += _b -> mass;
    velocity += _b -> velocity;
    temp = (temp + _b -> temp) / 2;
}

void Node3d::add_twob(Body3d* _b, Node3d* _n){
    assert(this -> is_empty());
    mass = _b -> mass + _n -> mass;
    velocity = _b -> velocity + _n -> velocity;
    temp = (_b -> temp + _n -> temp)/2;
    center_of_mass = (_b -> position * _b -> mass  + _n -> center_of_mass * _n -> mass) / mass;
}

/*                   //===---------------------------------------------------------===//
.                    // BARNES HUT IMPLEMENTATION                                     //
.                    //===---------------------------------------------------------===*/

Barnes_Hut3d::Barnes_Hut3d(sim3d_meta m){
    this -> data = m;
    data.type = BARNES_HUT;
    get_force = force_funcs_3d[data.selected_force];
    update_body = update_funcs_3d[data.selected_update];
    data.graph_height *= data.HEIGHT;

    get_rndX = std::uniform_real_distribution<double>(0, data.HEIGHT);
    get_rndY = std::uniform_real_distribution<double>(0, data.WIDTH);

    threads.reserve(data.NUM_THREADS);
    obj = bodies.size(); 
}
Barnes_Hut3d::Barnes_Hut3d(){
    this -> data = sim3d_meta();
    data.type = BARNES_HUT;
    get_force = force_funcs_3d[data.selected_force];
    update_body = update_funcs_3d[data.selected_update];
    data.graph_height *= data.HEIGHT;

    get_rndX = std::uniform_real_distribution<double>(0, data.HEIGHT);
    get_rndY = std::uniform_real_distribution<double>(0, data.WIDTH);

    threads.reserve(data.NUM_THREADS);
    obj = bodies.size(); 
}



void Barnes_Hut3d::step(){
    make_sections();
    make_tree();
    calculate_com();        

    update_bodies();
    nodes3d.clear();
    //std::cout << "stepped!";
}

//===---------------------------------------------------------===//
// Inserting body into tree                                      //
//===---------------------------------------------------------===//


bool Barnes_Hut3d::is_inbound(vec3& v){
    return 
        v.x < this -> data.WIDTH && 
        v.y < this -> data.WIDTH && 
        v.z < this -> data.WIDTH
    ;
}

void Barnes_Hut3d::make_sections(){
    sections.clear();
    int mult = bodies.size() / data.NUM_THREADS;

    for (int i = 0; i < data.NUM_THREADS-1; ++i)
        sections.emplace_back(i*mult);

    sections.emplace_back(bodies.size());
}

void Barnes_Hut3d::insert(Body3d* body){
    if (!is_inbound(body -> position)) return;

    int node_index = get_to_best_leaf(body);
    size_t cnode = node_index;

    if (nodes3d[cnode].is_empty()){
        nodes3d[cnode].init(body); //* it gets initialized
        return;
    }
    //* from now on we know it is not empty therefore it has a body ptr != nullptr

    if ((body -> position - nodes3d[cnode].center_of_mass).sqr_magn() < 100){ // are they in the same place?
        nodes3d[cnode].merge(body);
        return;
    }

    int tries = 5, new_node = node_index;
    size_t child = 0;
    short opt_1, opt_2;

    while (tries--){// trova un modo migliore
        // subdivides the current node
        child = subdivide(node_index);

        // gets where each body (the already existing one and the one we want to insert) wants to be
        opt_1 = opt_position(nodes3d[node_index].center_of_mass, nodes3d[node_index].center);
        opt_2 = opt_position(body -> position, nodes3d[node_index].center);

        if (opt_1 == opt_2) { // conflict
            new_node = child + opt_1;
            continue;
        }

        size_t
            p1 = child + opt_1,
            p2 = child + opt_2
        ;

        nodes3d[p1].init(&nodes3d[cnode]);
        nodes3d[p2].init(body);
        return;
    }

    nodes3d[child].add_twob(body, &nodes3d[cnode]);

}

int Barnes_Hut3d::get_to_best_leaf(Body3d* _b){
    /*
    * NOTE: there HAS to be a leaf node at the end of the vector
    * because every time a node is created it is initialized 
    * with "is_leaf = true" and if it's a leaf it does not have childrens
    * halting the loop.
    */
    size_t node_index = 0;

    // make sure the current node is a leaf
    while (!nodes3d[node_index].is_leaf()) {
        short displacement = opt_position(_b -> position, nodes3d[node_index].center);
        node_index = nodes3d[node_index].child + displacement; // skips to the next optimal node
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
int Barnes_Hut3d::opt_position(vec3 p, vec3 center){ 
    bool 
        x = p.x > center.x,
        y = p.y > center.y,
        z = p.z > center.z
    ;

    return x + y*2 + z*4;
}

   
void Barnes_Hut3d::calculate_com(){
    for (int i = nodes3d.size()-1; i >= 0; i--){
        if (nodes3d[i].is_leaf()) continue;

        nodes3d[i].mass = 0;
        nodes3d[i].center_of_mass = {0,0,0};

        for (int dis = 0; dis < 8; ++dis){
            size_t index = nodes3d[i].child +dis;

            auto& child = nodes3d[index];

            nodes3d[i].mass += child.mass;
            nodes3d[i].center_of_mass += child.center_of_mass*child.mass;
        }

        nodes3d[i].center_of_mass /= nodes3d[i].mass;
    }
}

void Barnes_Hut3d::get_node_body(size_t node, Body3d* body){
    assert(node < nodes3d.size());
    auto& cnode = nodes3d[node];
    
    if (cnode.is_empty()) return;

    double d_squared =(body -> position - nodes3d[node].center_of_mass).sqr_magn();

    if (d_squared * theta * theta > nodes3d[node].size * nodes3d[node].size){ // use the optmisation

        body -> acceleration += get_force(
            body -> mass, nodes3d[node].mass,
            body -> velocity, nodes3d[node].velocity,
            body -> position, nodes3d[node].center_of_mass,
            this
        ) / body -> mass;

        return;
    } 
    
    if (cnode.is_leaf()){

        if (nodes3d[node].center_of_mass == body -> position || !cnode.mass) return;

        body -> acceleration += get_force(
            body -> mass, nodes3d[node].mass,
            body -> velocity, nodes3d[node].velocity,
            body -> position, nodes3d[node].center_of_mass,
            this
        ) / body -> mass;

            
        return;
    }
  
    for (int i = 0; i < 8; i++)
        get_node_body(nodes3d[node].child + i, body);
    
}

void Barnes_Hut3d::update_bodies(){
    this -> obj = bodies.size(); 
    for (int i = 0; i < this -> data.NUM_THREADS; ++i)
        this -> threads.emplace_back(update_bundle, this, i);

    for (auto& t : threads)
        t.join();

    this -> threads.clear();
}

void Barnes_Hut3d::update_bundle(Barnes_Hut3d* _s, unsigned short index){
    unsigned int mult, start, stop;
    mult = _s -> obj/_s -> data.NUM_THREADS;
    start = index * mult;
    stop = (index + 1) * mult;
              
    for (int i = start; i < stop; ++i){  
        Body3d* body = &_s -> bodies[i];
        _s -> get_node_body(0, body);
        _s -> update_body(body, _s);

        body -> temp -= body -> acceleration.sqr_magn() ;//+ body -> temp * body -> temp * body -> temp * body -> temp *  _s -> data.dt;
    }
}

double variation(){
    return ((std::rand() % 100) + 50)/150; 
}

vec3 Barnes_Hut3d::pick_newpos(){
    return {
        data.WIDTH/2 * variation(),
        data.HEIGHT/2 * variation(),
        data.depth/2 *variation()
    };
}

void Barnes_Hut3d::make_tree(){
    //nodes3d.reserve(obj * 3);
    threads.clear();

    nodes3d.emplace_back(pick_newpos(), data.WIDTH);
    
    
    for (auto& b : bodies){
        insert(&b);
    }
}

size_t Barnes_Hut3d::subdivide(int n){
    //std::cout << nodes3d.size() << "\n";
    size_t child_index = nodes3d.size();
    Node3d& cnode = nodes3d[n]; 

    cnode.child = child_index;

    for (int i = 0; i < 8; i++){
       nodes3d.emplace_back(
           vec3({
           /*x = */ (0.5 - (i % 2 == 0)) * cnode.size + cnode.center.x,
           /*y = */ (0.5 - (i < 2 || (i > 3 && i < 5))) * cnode.size + cnode.center.y,
           /*z = */ (0.5 - (i < 3)) * cnode.size + cnode.center.z 
           }),
           cnode.size/2
       );
    }

    return child_index;
}



}
}