#pragma once

#include <vector>
#include <cassert>
#include <thread>
#include <type_traits>
#include <algorithm>
#include <bitset>
#include <unordered_map>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>

using namespace tbb;


#include "Aster/simulations/sim_obj.h"
#include "Aster/building-api/builder.h"

#include "Aster/graphics/2d_graphics.h"

#include "Aster/physics/vectors.h"
#include "Aster/physics/body.h" 

#include "Aster/simulations/barnes-hut.h"

#include <bitset>

namespace Aster{   

namespace GPU{
    void sort(uint64_t* input, uint64_t* output, size_t N);
}

template <typename T> void get_new_temp(Simulation<T>* _s, size_t body, T pos, T vel, REAL temp, REAL mass);

namespace Barnes{

template <typename F>
FORCE_INLINE void parallel(int cores, size_t num, F func) {

    const size_t n_threads = std::min(num, static_cast<size_t>(cores));

    if (num == 0) return;

    std::vector<std::thread> threads;
    threads.reserve(n_threads);

    const size_t chunk_size = (num + n_threads - 1) / n_threads;

    for (size_t i = 0; i < n_threads; ++i) {
        const size_t start = i * chunk_size;
        const size_t end = std::min(start + chunk_size, num);

        if (start >= end) break;

        threads.emplace_back([start, end, &func]() {
            for (size_t b = start; b < end; ++b) {
                func(b);
            }
        });
    }

    for (auto& t : threads) {
        t.join();
    }
}



/*                   //===---------------------------------------------------------===//
.                    // NODE IMPLEMENTATION                                           //
.                    //===---------------------------------------------------------===*/

template <typename T>
void NodesArray<T>::clear(){
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
template <typename T>
bool NodesArray<T>::is_empty(size_t i) const {
    return !masses[i] && is_leaf(i);
}

/**
* @brief returns if the node is empty
* @returns mass == 0 && is_leaf()
*/
template <typename T>
bool NodesArray<T>::is_leaf(size_t i) const {
    return right_nodes[i] == -1 && left_nodes[i] == -1;
}

 /**
 * @brief allocates the right amount of memory to the nodes
 */
template <typename T>
void NodesArray<T>::allocate(size_t n){
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
template <typename T>
void NodesArray<T>::resize(size_t n){
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

template <typename T>
size_t NodesArray<T>::size() const{
    return masses.size();
}

/**
* @brief initializes the node with the given body
* @param _b: ptr to the body to copy
*/
template <typename T>
void NodesArray<T>::init(Simulation<T>* _s,size_t  _b, size_t index){
    masses[index]         = _s -> bodies.get_mass_of(_b);
    centers_of_mass[index] = _s -> bodies.get_position_of(_b);
    velocities[index]      = _s -> bodies.get_velocity_of(_b);
    temps[index]           = _s -> bodies.get_temp_of(_b);
    left_nodes[index] = -1;
    right_nodes[index] = -1;
}

/**
* @brief pushes back a node with the given body
* @param _b: ptr to the body to copy
*/
template <typename T>
void NodesArray<T>::push_back(Simulation<T>* _s,size_t  _b){
    masses.push_back(_s -> bodies.get_mass_of(_b));
    centers_of_mass.push_back(_s -> bodies.get_position_of(_b));
    velocities.push_back(_s -> bodies.get_velocity_of(_b));
    temps.push_back(_s -> bodies.get_temp_of(_b));
    right_nodes.push_back(-1);
    left_nodes.push_back(-1);
}


/**
* @brief merges two nodes in one node
*/
template <typename T>
void NodesArray<T>::unite(size_t start){
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
template <typename T>
void NodesArray<T>::merge(Simulation<T>* _s, size_t _b, size_t index){
    // combines the key values
    REAL body_mass = _s -> bodies.get_mass_of(_b);
    REAL total_mass = masses[index] + body_mass;
    
    if (total_mass <= 0) return;

    centers_of_mass[index] = (_s -> bodies.get_position_of(_b) * body_mass + centers_of_mass[index] * masses[index]) / total_mass; 
    velocities[index] = (velocities[index] * masses[index] +  _s -> bodies.get_velocity_of(_b) * body_mass) / total_mass;
    masses[index] = total_mass;
    
    // for the temperature it uses the maximum
    temps[index] = std::max(REAL(temps[index]), _s -> bodies.get_temp_of(_b));
}

template <typename T>
void barnes_update_forces(Simulation<T>* _sim){
    
    auto _s = reinterpret_cast<Barnes_Hut<T>*>(_sim);
    _s ->make_sections(); // for threading
    _s -> make_tree(); // generates the quad-tree / oct-tree
        
    // calculates the forces
    parallel_for(size_t(0), _s -> bodies.positions.size(), [&, _s](size_t i){
        _s -> bodies.get_acc_of(i).reset();
        _s -> get_node_body(_s -> compressed_mortons_size, i, _s -> bounding_box.magnitude());
    });


}


/*                   //===---------------------------------------------------------===//
.                    // BARNES HUT IMPLEMENTATION                                     //
.                    //===---------------------------------------------------------===*/



template <typename T>
Barnes_Hut<T>::Barnes_Hut(){
    static_assert(std::is_same<T, vec2>::value || std::is_same<T, vec3>::value, "Invalid type for class construction");

    num_childs = 4; // sets up the child number

    if constexpr (std::is_same<T, vec3>::value) // if its 3d it sets  it up for an oct-tree
        num_childs = 8;
    
    // sets up the simulation
    this -> data = sim_meta(); 
    this -> data.type = BARNES_HUT;
    this -> get_force = get_force_func<T>(NEWTON);
    this -> update_bodies = get_update_func<T>(this -> data.selected_update, this -> uses_GPU());
    this -> data.graph_height *= this -> data.size.y;
    this -> update_forces =  static_cast<void(*)(Simulation<T>*)>(barnes_update_forces<T>);

    this -> threads.reserve(this -> get_cores());
    this -> obj = this -> bodies.positions.size(); 
}

/**
* @brief steps the simulation forward
*/
template <typename T>
void Barnes_Hut<T>::step(){
    barnes_update_forces(this);
    this -> update_bodies(this);

    // triggers graph and steps the time
    this -> trigger_all_graphs();
    this -> time_passed++;
}

/**
* @brief generates an array for the bodies slicing
*/
template <typename T> 
void Barnes_Hut<T>::make_sections(){
    sections.clear();
    int mult = this -> bodies.positions.size() / this -> get_cores();

    for (int i = 0; i < this -> get_cores()-1; ++i)
        sections.emplace_back(i*mult);

    sections.emplace_back(this -> bodies.positions.size());
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
template <typename T>
void Barnes_Hut<T>::get_node_body(signed long int node, size_t index, REAL size){
    if (node < 0 || node >= static_cast<signed long int>(nodes.size()))
        return;

    REAL d_squared = (this -> bodies.get_position_of(index) - nodes.centers_of_mass[node]).sqr_magn();
    
    if (d_squared == 0) return; 

    if ( size * size < d_squared *theta * theta || nodes.is_leaf(node)){ // use the optmisation
        this -> bodies.get_acc_of(index) += this -> get_force(
            this -> bodies.get_mass_of(index), nodes.masses[node], 
            this -> bodies.get_velocity_of(index), nodes.velocities[node], 
            this -> bodies.get_position_of(index), nodes.centers_of_mass[node], 
            this
        ) / this -> bodies.get_mass_of(index);

        return;
    } 

    get_node_body(nodes.left_nodes[node],  index, size / 2);
    get_node_body(nodes.right_nodes[node], index, size / 2);
}

/**
* @brief internal function to update a batch of bodies
* @param _s: ptr to the simulation
* @param index: what chunck of the bodies array to scan
*/
template <typename T>
void update_bundle(Barnes_Hut<T>* _s, unsigned short index){
    if (_s -> nodes.size() <= 0) return;
    size_t N = _s->obj;
    size_t P = _s->get_cores();
    size_t chunk = (N + P - 1) / P;
    size_t start = index * chunk;
    size_t stop  = std::min(start + chunk, N);
    
    // calculates the forces of those bodies
    for (size_t i = start; i < stop; ++i){
        _s -> bodies.get_acc_of(i).reset();
        
        _s -> get_node_body(_s -> compressed_mortons_size, i, _s -> bounding_box.magnitude());
    }
}

template <typename T> 
void compute_tidal_heating(Barnes_Hut<T>* _s, size_t index){
    for (int i = _s -> nodes[0].child; i < _s -> num_childs + _s -> nodes[0].child; ++i){
        if (i > _s -> bodies.positions.size()) return;
        
        get_new_temp<T>(_s,
            index,_s ->  nodes.centers_of_mass[index], 
            _s -> nodes.velocities[index], _s -> nodes.temps[index], 
            _s -> nodes.masses[index]
        );
    }
}

inline signed int count_leading_zeros(uint32_t num){
    return num ? __builtin_clz(num) : 32;
};

template <typename T>
uint32_t get_tbreaked_morton(Barnes_Hut<T>* _s, const T& pos, size_t body_index) {
    return get_morton<T>(_s, pos);
}

template <typename T>
inline void translate2nodes(Simulation<T>* _s, NodesArray<T>& base_layer, const std::vector<std::pair<uint32_t, unsigned int>>& mortons){
    base_layer.clear();
    base_layer.resize(mortons.size());

    base_layer.init(_s,mortons[0].second, 0);
    size_t last = 0;

    for (size_t i = 1; i < mortons.size(); ++i) {
        if (mortons[i-1].first == mortons[i].first)
            base_layer.merge(_s, mortons[i].second, last);
        else{
            last++;
            base_layer.init(_s, mortons[i].second, last);
        }
    }
    base_layer.resize(last+1);
}

template <typename T>
Barnes_Hut<T>* Barnes_Hut<T>::set_theta(REAL _t){
    assert(_t > 0);
    theta = _t;
    return this;
}

template <typename T>
Simulation<T>* Barnes_Hut<T>::use_GPU(){
    this -> GPU_on = true;
    return this;
}

template <typename T>
void printTreeRecursive(const NodesArray<T>& nodes, int node, std::string prefix = "", bool isRight = false) {
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
template <typename T>
void print_nodes(const NodesArray<T>& nodes, int root = 0) {
    std::cout << "Tree structure:\n";
    printTreeRecursive(nodes, root);
}

template <typename T>
void Barnes_Hut<T>::make_tree(){
    this->threads.clear();
    this->nodes.clear();
    this->mortons.clear();

    size_t N = this->bodies.positions.size();
    if (N == 0) return;
    
    std::vector<std::pair<uint32_t, unsigned int>> enhanced_mortons;
    enhanced_mortons.reserve(N);
    this->bounding_box = this->get_center() * 2;

    for (size_t i = 0; i < N; ++i) {
        uint32_t morton_code = get_tbreaked_morton<T>(this, this->bodies.positions[i], i);
        enhanced_mortons.emplace_back(morton_code, i);
    }

    std::sort(enhanced_mortons.begin(), enhanced_mortons.end(), [](const std::pair<uint32_t, size_t>& a, const std::pair<uint32_t, size_t>&b ){ return a.first < b.first;});

    translate2nodes<T>(this, nodes, enhanced_mortons);
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

}}

}