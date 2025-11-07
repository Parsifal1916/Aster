
#include <vector>
#include <cassert>
#include <thread>
#include <type_traits>
#include <algorithm>
#include <bitset>
#include <unordered_map>
#include <bitset>
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
    #include <xmmintrin.h>
#endif

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

void NodesArray::resize(size_t n) {
    if (n == left_nodes.size()) return;
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

void Barnes_Hut::get_node_body(long root, size_t index, REAL size) {
    struct Task { long node; REAL size; };
    std::array<Task, 64> stack;
    int sp = 0;
    stack[sp++] = { root, size };

    auto& bodies = this->_s->bodies;
    const vec3& pos_i = bodies.positions[index];
    const vec3& vel_i = bodies.velocities[index];
    const REAL mass_i = bodies.masses[index];
    const REAL inv_m = REAL(1.0) / mass_i;
    const REAL theta2 = this->_s->theta * this->_s->theta;
    constexpr REAL half = REAL(0.5);
    vec3 n_acc;

    while (sp > 0) {
        const auto [node, s] = stack[--sp];
        if ((unsigned long)node >= nodes.size()) continue;

#if (defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86))
        _mm_prefetch(reinterpret_cast<const char*>(&nodes.centers_of_mass[node]), _MM_HINT_T0);
        _mm_prefetch(reinterpret_cast<const char*>(&nodes.velocities[node]), _MM_HINT_T0);
        _mm_prefetch(reinterpret_cast<const char*>(&nodes.masses[node]), _MM_HINT_T0);
        _mm_prefetch(reinterpret_cast<const char*>(&nodes.left_nodes[node]), _MM_HINT_T0);
        _mm_prefetch(reinterpret_cast<const char*>(&nodes.right_nodes[node]), _MM_HINT_T0);
#endif

        const vec3& com = nodes.centers_of_mass[node];
        const vec3& vel_node = nodes.velocities[node];
        const REAL mass_node = nodes.masses[node];

        const auto d = pos_i - com;
        const REAL d2 = d.sqr_magn();
        if (d2 == REAL(0)) continue;

        if ((s * s < d2 * theta2) || nodes.is_leaf(node)) {
            const auto acc = this->get_force(mass_i, mass_node, vel_i, vel_node, pos_i, com, this->_s);
            n_acc += acc * inv_m;
        } else {
            const long left  = nodes.left_nodes[node];
            const long right = nodes.right_nodes[node];

            if ((unsigned long)left < nodes.size()) {
#if (defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86))
                _mm_prefetch(reinterpret_cast<const char*>(&nodes.centers_of_mass[left]), _MM_HINT_T0);
#endif
                stack[sp++] = { left, s * half };

            }
            if ((unsigned long)right < nodes.size()) {
#if (defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86))
                _mm_prefetch(reinterpret_cast<const char*>(&nodes.centers_of_mass[right]), _MM_HINT_T0);
#endif
                stack[sp++] = { right, s * half };
            }
        }
    }

    bodies.accs[index] = n_acc;
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

    std::cout << prefix;
    if (!prefix.empty()) {
        std::cout << (isRight ? "└── " : "├── ");
    }
    std::cout << node << "\n";

    std::string newPrefix = prefix + (isRight ? "    " : "│   ");

    printTreeRecursive(nodes, nodes.left_nodes[node], newPrefix, false);
    printTreeRecursive(nodes, nodes.right_nodes[node], newPrefix, true);
}


void print_nodes(const NodesArray& nodes, int root = 0) {
    std::cout << "Tree structure:\n";
    printTreeRecursive(nodes, root);
}

void Barnes_Hut::make_tree(){
    std::vector<std::pair<uint32_t, unsigned int>> enhanced_mortons(get_range());
    
    size_t num_leaves = get_range();
    size_t num_internal = num_leaves - 1;
    if (num_leaves <= 1) return;
  
    this -> nodes.resize(num_internal+num_leaves);
    this->bounding_box = this->_s -> get_center() * 2;

    parallel_for(size_t(get_lower_bound()), size_t(get_upper_bound()), [&](size_t i){
        uint32_t morton_code = get_tbreaked_morton(this, this->_s->bodies.positions[i], i);
        enhanced_mortons[i - get_lower_bound()] = {morton_code, i - get_lower_bound()};
    });

    parallel_sort(enhanced_mortons.begin(), enhanced_mortons.end(), [](const std::pair<uint32_t, size_t>& a, const std::pair<uint32_t, size_t>&b ){ return a.first < b.first;});

    parallel_for(size_t(0), enhanced_mortons.size(), [&](size_t i){
        int _b = enhanced_mortons[i].second;
        auto& bodies = _s -> bodies;

        this->nodes.masses[i]         =  bodies.get_mass_of(_b);
        this->nodes.centers_of_mass[i] = bodies.get_position_of(_b);
        this->nodes.velocities[i]      = bodies.get_velocity_of(_b);
        this->nodes.temps[i]           = bodies.get_temp_of(_b);
        this->nodes.left_nodes[i] = -1;
        this->nodes.right_nodes[i] = -1;
    });

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
    size_t N = this -> _s -> bodies.positions.size();
    const size_t root = N;
    const REAL bd_size = (this -> _s -> get_center()*2).magnitude();
    this -> make_tree();

    N = this->get_range();
    
    if (warn_if(N <= 0, "either upper or lower bound were set too high in the gravity solver, resulting in a negative amount of bodies to load")) 
        return; 
        

    parallel_for(size_t(get_lower_bound()), size_t(get_upper_bound()), [&, this, root](size_t idx){
        this -> _s -> bodies.accs[idx].reset();
        this -> get_node_body(root, idx, bd_size);
    });
} 

}}