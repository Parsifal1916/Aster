#pragma once

#include <vector>
#include <cassert>
#include <thread>
#include <type_traits>
#include <algorithm>
#include <unordered_map>
#include <iterator>
#include <future>

#define CL_TARGET_OPENCL_VERSION 300
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif

#include "Aster/simulations/sim_obj.h"
#include "Aster/building-api/builder.h"

#include "Aster/graphics/2d_graphics.h"

#include "Aster/physics/vectors.h"
#include "Aster/physics/body.h"

#include "Aster/simulations/barnes_hut_GPU.h"
#include "Aster/building-api/GPU_endpoint.h"
#include "Aster/impl/kernels/barnes_tree_creation.cl.h"
#include "Aster/impl/kernels/bitonic_sort.cl.h"
#include "Aster/impl/kernels/basic_barnes_force.cl.h"

#include <bitset>

namespace Aster{
namespace Barnes{

#define check(call)                                                    \
do {                                                                   \
    cl_int _err = (call);                                              \
    if (_err != CL_SUCCESS) {                                          \
        fprintf(stderr, "OpenCL error %d at %s:%d – %s\n",             \
                _err, __FILE__, __LINE__, #call);                      \
        abort();                                                       \
    }                                                                  \
} while (0)

FORCE_INLINE void load_const_buff(size_t size, void* data, cl_mem& buff){
    using namespace GPU;
    cl_int operation_result;

    buff = clCreateBuffer(context, CL_MEM_READ_ONLY, size, nullptr, &operation_result);
    check(operation_result);
    operation_result = clEnqueueWriteBuffer(queue, buff, CL_FALSE, 0, size, data, 0, nullptr, nullptr);
    check(operation_result);
}

template <typename T> 
void BHG<T>::upload_force_calc(size_t num_leaves, size_t tree_size, size_t num_bodies){
    using namespace GPU;
    cl_int operation_result = 0 ;
    if (num_leaves == 0) return;
    int mult = sizeof(T) / sizeof(REAL);
    
    size_t double2_size_tree = sizeof(REAL) * mult * tree_size;
    size_t int_size_tree = sizeof(int) * tree_size;
    size_t double_size_tree = sizeof(REAL) * tree_size;

    size_t double2_size_bodies = sizeof(REAL) * mult * num_bodies; 
    size_t double_size_bodies = sizeof(REAL) * num_bodies;       
    
    cl_mem com_buff, rights_buff, lefts_buff, nmasses_buff, pos_buff, bmasses_buff;

    load_const_buff(double2_size_tree,   this->nodes.centers_of_mass.data(), com_buff);       
    load_const_buff(int_size_tree,       this->nodes.left_nodes.data(),      lefts_buff);     
    load_const_buff(int_size_tree,       this->nodes.right_nodes.data(),     rights_buff);    
    load_const_buff(double_size_tree,    this->nodes.masses.data(),          nmasses_buff);   
    
    load_const_buff(double2_size_bodies, this->bodies.positions.data(),      pos_buff);       
    load_const_buff(double_size_bodies,  this->bodies.masses.data(),         bmasses_buff);   


    // accelerations
    cl_mem accs_buff = clCreateBuffer(context,
            CL_MEM_WRITE_ONLY,
            double2_size_bodies,
            nullptr,
        &operation_result);
    check(operation_result);

    size_t LW_size = 64;
    size_t GW_size = ((num_bodies + LW_size - 1) / LW_size) * LW_size;

    int root_node = num_leaves;  
    REAL tree_bound_size = this->bounding_box.magnitude(); 
    
    operation_result = clSetKernelArg(force_calculator, 0, sizeof(int), &num_bodies);
    check(operation_result);
    operation_result = clSetKernelArg(force_calculator, 1, sizeof(REAL), &this->data.G);
    check(operation_result);
    operation_result = clSetKernelArg(force_calculator, 2, sizeof(REAL), &this->theta);
    check(operation_result);
    operation_result = clSetKernelArg(force_calculator, 3, sizeof(REAL), &tree_bound_size);  
    check(operation_result);
    operation_result = clSetKernelArg(force_calculator, 4, sizeof(int), &root_node);
    check(operation_result);
    operation_result = clSetKernelArg(force_calculator, 5, sizeof(int), &tree_size); 
    check(operation_result);
    operation_result = clSetKernelArg(force_calculator, 6, sizeof(cl_mem), &nmasses_buff);
    check(operation_result);
    operation_result = clSetKernelArg(force_calculator, 7, sizeof(cl_mem), &pos_buff); 
    check(operation_result);
    operation_result = clSetKernelArg(force_calculator, 8, sizeof(cl_mem), &com_buff); 
    check(operation_result);
    operation_result = clSetKernelArg(force_calculator, 9, sizeof(cl_mem), &lefts_buff); 
    check(operation_result);
    operation_result = clSetKernelArg(force_calculator, 10, sizeof(cl_mem), &rights_buff); 
    check(operation_result);
    operation_result = clSetKernelArg(force_calculator, 11, sizeof(cl_mem), &bmasses_buff);  
    check(operation_result);
    operation_result = clSetKernelArg(force_calculator, 12, sizeof(cl_mem), &accs_buff);  
    check(operation_result);


    operation_result = clEnqueueNDRangeKernel(queue, force_calculator, 1, 0, &GW_size, &LW_size, 0, nullptr, nullptr);
    check(operation_result);


    operation_result = clEnqueueReadBuffer(queue, accs_buff,
        CL_FALSE, 0, 
        double2_size_bodies,
        this->bodies.accs.data(),
        0, nullptr, nullptr
    ); 
    check(operation_result); 

    clFinish(queue);
    clReleaseMemObject(com_buff);
    clReleaseMemObject(rights_buff);
    clReleaseMemObject(lefts_buff); 
    clReleaseMemObject(nmasses_buff); 
    clReleaseMemObject(pos_buff); 
    clReleaseMemObject(bmasses_buff);
    clReleaseMemObject(accs_buff);

}

template <typename T> 
void BHG<T>::load_tree_kernel(int n, int num_leaves, void* mortons){
    using namespace GPU;
    cl_int operation_result;
    if (num_leaves == 0) return;
    
    size_t children_size = (n + num_leaves) * sizeof(int);
    size_t mortons_size = num_leaves * sizeof(cl_uint2); 
 
    cl_mem left_nodes_buff  = clCreateBuffer(context, CL_MEM_READ_WRITE, children_size, nullptr, &operation_result);
    check(operation_result);
    cl_mem right_nodes_buff = clCreateBuffer(context, CL_MEM_READ_WRITE, children_size, nullptr, &operation_result);
    check(operation_result);
    cl_mem mortons_buff     = clCreateBuffer(context, CL_MEM_READ_ONLY,  mortons_size,  nullptr, &operation_result);
    check(operation_result);
    
    check(clEnqueueWriteBuffer(queue, mortons_buff, CL_FALSE, 0, mortons_size, 
                                  mortons, 0, nullptr, nullptr));
    check(clEnqueueWriteBuffer(queue, left_nodes_buff,  CL_FALSE, 0, children_size,
                                  this->nodes.left_nodes.data(), 0, nullptr, nullptr));
    check(clEnqueueWriteBuffer(queue, right_nodes_buff, CL_FALSE, 0, children_size,
                                  this->nodes.right_nodes.data(), 0, nullptr, nullptr));


    size_t LW_size = 64;
    size_t GW_size = ((n + LW_size - 1) / LW_size) * LW_size;

    operation_result = clSetKernelArg(tree_builder, 0, sizeof(int), &num_leaves);
    check(operation_result);
    operation_result = clSetKernelArg(tree_builder, 1, sizeof(cl_mem), &mortons_buff);
    check(operation_result);
    operation_result = clSetKernelArg(tree_builder, 2, sizeof(cl_mem), &left_nodes_buff);
    check(operation_result);
    operation_result = clSetKernelArg(tree_builder, 3, sizeof(cl_mem), &right_nodes_buff);
    check(operation_result);


    operation_result = clEnqueueNDRangeKernel(queue, tree_builder, 1, 0, &GW_size, &LW_size, 0, nullptr, nullptr );
    check(operation_result);

    cl_event evt[2];
    check(clEnqueueReadBuffer(queue, right_nodes_buff, CL_FALSE, 0, children_size,
                                 this->nodes.right_nodes.data(), 0, nullptr, &evt[0]));
    check(clEnqueueReadBuffer(queue, left_nodes_buff,  CL_FALSE, 0, children_size,
                                 this->nodes.left_nodes.data(),  0, nullptr, &evt[1]));
    
    check(clWaitForEvents(2, evt));
    clReleaseEvent(evt[0]); clReleaseEvent(evt[1]);
    
    check(clFinish(queue));
    clReleaseMemObject(left_nodes_buff);
    clReleaseMemObject(right_nodes_buff);
    clReleaseMemObject(mortons_buff);

}

template <typename T>
inline void translate2nodes(Simulation<T>* _s, NodesArray<T>& base_layer, const std::vector<cl_uint2>& mortons){
    base_layer.clear();
    base_layer.resize(mortons.size());

    base_layer.init(_s,mortons[0].s[1], 0);
    size_t last = 0;

    for (size_t i = 1; i < mortons.size(); ++i) {
        if (mortons[i-1].s[0] == mortons[i].s[0])
            base_layer.merge(_s, mortons[i].s[1], last);
        else{
            last++;
            base_layer.init(_s, mortons[i].s[1], last);
        }
    }

    base_layer.resize(last+1);
}

template <typename T> 
void BHG<T>::build_tree(){
    this->threads.clear();
    this->nodes.clear();
    this->mortons.clear();

    size_t N = this->bodies.positions.size();
    if (N == 0) return;
    
    std::vector<std::pair<uint32_t, unsigned int>> enhanced_mortons;
    enhanced_mortons.resize(N);
    this->bounding_box = this->get_center() * 2;

    parallel_for(size_t(0), N, [this, &enhanced_mortons](size_t i){
        uint32_t morton_code = get_tbreaked_morton<T>(this, this->bodies.positions[i], i);
        enhanced_mortons[i] = {morton_code, i};
    });

    parallel_sort(enhanced_mortons.begin(), enhanced_mortons.end(), [](const std::pair<uint32_t, size_t>& a, const std::pair<uint32_t, size_t>&b ){ return a.first < b.first;});

    translate2nodes<T>(this, this ->  nodes, enhanced_mortons);
    size_t num_leaves = this -> nodes.size();
    
    
    if (num_leaves <= 1) return;
  
    size_t num_internal = num_leaves - 1;
    this -> nodes.resize(num_leaves + num_internal);
    
    this -> load_tree_kernel(num_internal, num_leaves, enhanced_mortons.data());
    
    this -> nodes.unite(num_leaves);
    
    //for (int i = 0; i < this -> nodes.size(); ++i){
    //    std::cout << "Node with mass  : " << this -> nodes.masses[i] << "\n"
    //              << "Node with center: ("<< this -> nodes.centers_of_mass[i].x << ", " << this -> nodes.centers_of_mass[i].y << ")\n"
    //              << "Left            : " << this -> nodes.left_nodes[i] << "\n"
    //              << "Right           : " << this -> nodes.right_nodes [i] << "\n"
    //              << "==============================\n";
    //}

    this->compressed_mortons_size = num_leaves;
}

template <typename T>
BHG<T>::BHG() {
    // sets up the simulation
    this -> data = sim_meta(); 
    this -> data.type = BH_cl;
    this -> get_force = get_force_func<T>(NEWTON);
    this -> update_bodies = get_update_func<T>(this -> data.selected_update, this -> uses_GPU());
    this -> data.graph_height *= this -> data.size.y;
    this -> update_forces =  static_cast<void(*)(Simulation<T>*)>(barnes_update_forces<T>);

    this -> threads.reserve(this -> get_cores());
    this -> obj = this -> bodies.positions.size(); 

    GPU::init_opencl();
    std::string tree_kernel_name = "build_tree";
    std::string force_kernel_name = "barnes_force";
    //std::string lbit_kernel_name = "local_block_bitonic";
    //std::string gbit_kernel_name = "global_bitonic_merge";
    GPU::compile_kernel(&tree_kernel_name, &GPU::barnes_tree_cl, tree_builder);
    
    //if constexpr (sizeof(REAL) == sizeof(double)) 
    
    if constexpr (std::is_same<T, vec2>::value)
        GPU::compile_kernel(&force_kernel_name, &GPU::barnes_force_basic, force_calculator);
    else
        GPU::compile_kernel(&force_kernel_name, &GPU::barnes_force_basic_3d, force_calculator);
        //else 
    //    GPU::compile_kernel(&force_kernel_name, &GPU::barnes_force_basic_lite, force_calculator);
    //GPU::compile_kernel(&lbit_kernel_name, &GPU::bitonic_sort_cl, this -> local_block_bitonic);
    //GPU::compile_kernel(&gbit_kernel_name, &GPU::bitonic_sort_cl, this -> global_bitonic_merge);
}


template <typename T>
void BHG<T>::step(){
    this -> make_sections(); // for threading
    this -> build_tree(); // generates the quad-tree / oct-tree
    
    upload_force_calc(this -> compressed_mortons_size, this ->nodes.size(), this -> bodies.positions.size());
    //barnes_update_forces(this);

    this -> update_bodies(this);

    // triggers graph and steps the time
    this -> trigger_all_graphs();
    this -> time_passed++;
}

}   
}