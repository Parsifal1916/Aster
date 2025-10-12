
#include <vector>
#include <cassert>
#include <thread>
#include <type_traits>
#include <algorithm>
#include <bitset>
#include <unordered_map>
#include <bitset>

#include "Aster/simulations/barnes-hut.h"
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>

using namespace tbb;




namespace Aster{   

 void get_new_temp(Simulation* _s, size_t body, vec3 pos, vec3 vel, REAL temp, REAL mass);

namespace Barnes{



/*                   //===---------------------------------------------------------===//
.                    // NODE IMPLEMENTATION                                           //
.                    //===---------------------------------------------------------===*/


void NodesArray::clear(){
    centers_of_mass.clear();
    temps.clear();
    velocities.clear();
    right_nodes.clear();
    left_nodes.clear();
    masses.clear();
}

/**
* @brief returns if the node is empty
* @returns mass == 0 && is_leaf()
*/

bool NodesArray::is_empty(size_t i) const {
    return !masses[i] && is_leaf(i);
}

/**
* @brief returns if the node is empty
* @returns mass == 0 && is_leaf()
*/

bool NodesArray::is_leaf(size_t i) const {
    return right_nodes[i] == -1 && left_nodes[i] == -1;
}

 /**
 * @brief allocates the right amount of memory to the nodes
 */

void NodesArray::allocate(size_t n){
    centers_of_mass.reserve(n);
    temps.reserve(n);
    velocities.reserve(n);
    right_nodes.reserve(n);
    left_nodes.reserve(n);
    masses.reserve(n);
}

/**
@brief resizes all arrays
*/

void NodesArray::resize(size_t n){
    size_t prev_size =  left_nodes.size();
    centers_of_mass.resize(n);
    temps.resize(n);
    velocities.resize(n);
    right_nodes.resize(n);
    left_nodes.resize(n);
    masses.resize(n);
    for (int i = prev_size; i < n; ++i){
        left_nodes[i] = -1;
        right_nodes[i] = -1;
    }
}


size_t NodesArray::size() const{
    return masses.size();
}

/**
* @brief initializes the node with the given body
* @param _b: ptr to the body to copy
*/

void NodesArray::init(Simulation* _s,size_t  _b, size_t index){
    auto& bodies = _s -> bodies;

    masses[index]         =  bodies.get_mass_of(_b);
    centers_of_mass[index] = bodies.get_position_of(_b);
    velocities[index]      = bodies.get_velocity_of(_b);
    temps[index]           = bodies.get_temp_of(_b);
    left_nodes[index] = -1;
    right_nodes[index] = -1;
}

/**
* @brief pushes back a node with the given body
* @param _b: ptr to the body to copy
*/

void NodesArray::push_back(Simulation* _s,size_t  _b){
    auto& bodies = _s -> bodies;
    
    masses.push_back(bodies.get_mass_of(_b));
    centers_of_mass.push_back(bodies.get_position_of(_b));
    velocities.push_back(bodies.get_velocity_of(_b));
    temps.push_back(bodies.get_temp_of(_b));
    right_nodes.push_back(-1);
    left_nodes.push_back(-1);
}


/**
* @brief merges two nodes in one node
*/

void NodesArray::unite(size_t start){
    if (left_nodes[start] == -1 || right_nodes[start] == -1) return;
    
    if (left_nodes[start] < 0 || left_nodes[start] >= size() || right_nodes[start] < 0 || right_nodes[start] >= size()) return;

    unite(left_nodes[start]);
    unite(right_nodes[start]);
    
    REAL total_mass = masses[right_nodes[start]] + masses[left_nodes[start]];
    if (total_mass <= 0) return;

    masses[start] = total_mass;
    centers_of_mass[start] = (centers_of_mass[right_nodes[start]] * masses[right_nodes[start]] + centers_of_mass[left_nodes[start]] * masses[left_nodes[start]]) / masses[start]; 
    velocities[start] = (velocities[right_nodes[start]] * masses[right_nodes[start]] + velocities[left_nodes[start]] * masses[left_nodes[start]]) / masses[start];
    temps[start] = std::max(temps[right_nodes[start]], temps[left_nodes[start]]);
}

/**
* @brief merges two bodies in one node in case of conflict
* @param body: index to the body to merge
*/

void NodesArray::merge(Simulation* _s, size_t _b, size_t index){
    // combines the key values
    auto& bodies = _s -> bodies;
    REAL body_mass = bodies.get_mass_of(_b);
    REAL total_mass = masses[index] + body_mass;
    
    if (total_mass <= 0) return;

    centers_of_mass[index] = (bodies.get_position_of(_b) * body_mass + centers_of_mass[index] * masses[index]) / total_mass; 
    velocities[index] = (velocities[index] * masses[index] +  bodies.get_velocity_of(_b) * body_mass) / total_mass;
    masses[index] = total_mass;
    
    // for the temperature it uses the maximum
    temps[index] = std::max(REAL(temps[index]), bodies.get_temp_of(_b));
}


/*                   //===---------------------------------------------------------===//
.                    // BARNES HUvec3 IMPLEMENTATION                                     //
.                    //===---------------------------------------------------------===*/




Barnes_Hut::Barnes_Hut(Simulation* _s){
    num_childs = 4; 
    
    this -> _s = _s;
}



//===---------------------------------------------------------===//
// Inserting body into tree                                      //
//===---------------------------------------------------------===//


/**
* @brief calculates the force acting between a node and a body **recursive**
* @param node: node to calculate the force from
* @param body: ptr to the body being acted upon
* @returns nothing, everything is done internally
*/

void Barnes_Hut::get_node_body(signed long int node, size_t index, REAL size){
    if (node < 0 || node >= static_cast<signed long int>(nodes.size()))
        return;

    auto& bodies = this -> _s -> bodies;

    REAL d_squared = (bodies.get_position_of(index) - nodes.centers_of_mass[node]).sqr_magn();
    
    if (d_squared == 0) return; 

    if ( size * size < d_squared *theta * theta || nodes.is_leaf(node)){ // use the optmisation
        bodies.get_acc_of(index) += this -> get_force(
            bodies.get_mass_of(index), nodes.masses[node], 
            bodies.get_velocity_of(index), nodes.velocities[node], 
            bodies.get_position_of(index), nodes.centers_of_mass[node], 
            this -> _s
        ) / bodies.get_mass_of(index);

        return;
    } 

    get_node_body(nodes.left_nodes[node],  index, size / 2);
    get_node_body(nodes.right_nodes[node], index, size / 2);
}

inline signed int count_leading_zeros(uint32_t num){
    return num ? __builtin_clz(num) : 32;
};


uint32_t get_tbreaked_morton(Barnes_Hut* _s, const vec3& pos, size_t body_index) {
    return get_morton(_s, pos);
}


inline void translate2nodes(Simulation* _s, NodesArray& base_layer, const std::vector<std::pair<uint32_t, unsigned int>>& mortons){
    base_layer.clear();
    base_layer.resize(mortons.size());

    base_layer.init(_s,mortons[0].second, 0);
    size_t last = 0;

    for (size_t i = 1; i < mortons.size(); ++i) {
        last++;
        base_layer.init(_s, mortons[i].second, last);
    }
    base_layer.resize(last+1);
}


void printTreeRecursive(const NodesArray& nodes, int node, std::string prefix = "", bool isRight = false) {
    if (node == -1) return;

    // Stampa il nodo corrente con connettori
    std::cout << prefix;
    if (!prefix.empty()) {
        std::cout << (isRight ? "└── " : "├── ");
    }
    std::cout << node << "\n";

    // Aggiorna prefisso per i figli
    std::string newPrefix = prefix + (isRight ? "    " : "│   ");

    // Stampa ramo sinistro e destro
    printTreeRecursive(nodes, nodes.left_nodes[node], newPrefix, false);
    printTreeRecursive(nodes, nodes.right_nodes[node], newPrefix, true);
}

// Wrapper

void print_nodes(const NodesArray& nodes, int root = 0) {
    std::cout << "Tree structure:\n";
    printTreeRecursive(nodes, root);
}


void Barnes_Hut::make_tree(){
    this->threads.clear();
    this->nodes.clear();
    this->mortons.clear();

    size_t N = this->_s -> bodies.positions.size();
    if (N == 0) return; 
    
    std::vector<std::pair<uint32_t, unsigned int>> enhanced_mortons;
    enhanced_mortons.reserve(N);
    this->bounding_box = this->_s -> get_center() * 2;

    for (size_t i = 0; i < N; ++i) {
        uint32_t morton_code = get_tbreaked_morton(this, this->_s -> bodies.positions[i], i);
        enhanced_mortons.emplace_back(morton_code, i);
    }

    std::sort(enhanced_mortons.begin(), enhanced_mortons.end(), [](const std::pair<uint32_t, size_t>& a, const std::pair<uint32_t, size_t>&b ){ return a.first < b.first;});

    translate2nodes(this -> _s, nodes, enhanced_mortons);
    size_t num_leaves = nodes.size();
    
    if (num_leaves <= 1) return;
  
    size_t num_internal = num_leaves - 1;
    nodes.resize(num_leaves + num_internal);

    auto delta = [&](size_t i, size_t j) -> int {
        if (j >= enhanced_mortons.size()) return -1;
        if (enhanced_mortons[i].first == enhanced_mortons[j].first) {
            return 32 + __builtin_clz(i ^ j); 
        }
        return __builtin_clz(enhanced_mortons[i].first ^ enhanced_mortons[j].first);
    };

    parallel_for(size_t(0), num_internal, [&, delta, this](size_t i){
        int d_left = (i == 0) ? -1 : delta(i, i - 1);
        int d_right = delta(i, i + 1);
        int d = (d_right > d_left) ? 1 : -1;
        
        int d_min = (d == 1) ? d_left : d_right;
        size_t l_max = 2;
        while ((int)(i + l_max * d) >= 0 && i + l_max * d < num_leaves && delta(i, i + l_max * d) > d_min) {
            l_max *= 2;
        }
        
        size_t l = 0;
        for (size_t t = l_max / 2; t >= 1; t /= 2) {
            size_t test_pos = i + (l + t) * d;
            if ((int)test_pos >= 0 && test_pos < num_leaves && delta(i, test_pos) > d_min) {
                l += t;
            }
        }
        
        size_t j = i + l * d;
        int d_node = delta(i, j);
        
        size_t s = 0;
        size_t range_size = (j > i) ? j - i : i - j;
        for (size_t t = (range_size + 1) / 2; t >= 1; t = (t + 1) / 2) {
            size_t test_pos = i + (s + t) * d;
            if ((int)test_pos >= 0 && test_pos < num_leaves && delta(i, test_pos) > d_node) {
                s += t;
            }
            if (t == 1) break;
        }
        
        size_t gamma = i + s * d + std::min(d, 0);
        
        size_t internal_node = num_leaves + i;
        
        size_t left_range = std::min(i, j);
        size_t right_range = std::max(i, j);
         
        if (left_range == gamma) {
            nodes.left_nodes[internal_node] = gamma; 
        } else {
            nodes.left_nodes[internal_node] = num_leaves + gamma; 
        }
        
        if (right_range == gamma + 1) {
            nodes.right_nodes[internal_node] = gamma + 1;  
        } else {
            nodes.right_nodes[internal_node] = num_leaves + gamma + 1; 
        }
    });


    nodes.unite(num_leaves);
    this->compressed_mortons_size = num_leaves;

}
void Barnes_Hut::compute_forces(){
    const size_t N = this -> _s -> bodies.positions.size();
    const REAL bd_size = (this -> _s -> get_center()*2).magnitude();
    this -> make_tree();

    parallel_for(size_t(0), N, [&, this, N](size_t idx){
        this -> get_node_body(N, idx, bd_size);
    });
}

}}