#include <vector>
#include <cassert>
#include <thread>
#include <type_traits>
#include <algorithm>
#include <unordered_map>
#include <iterator>
#include <future>

#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>
using namespace tbb;
#define CL_TARGET_OPENCL_VERSION 300
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif

#include "Aster/simulations/barnes_hut_GPU.h"

#include "Aster/building-api/GPU_endpoint.h"
#include "Aster/kernels/PNs.cl.h"
#include "Aster/kernels/PNs3d.cl.h"
#include "Aster/kernels/newton_update.cl.h"
#include "Aster/kernels/newton_update3d.cl.h"
#include "Aster/kernels/barnes_tree_creation.cl.h"
#include "Aster/kernels/basic_barnes_force.cl.h"

#include <bitset>

namespace Aster{
namespace Barnes{

FORCE_INLINE void load_const_buff(size_t size, void* data, cl_mem& buff){
    using namespace GPU;
    cl_int operation_result;

    buff = clCreateBuffer(context, CL_MEM_READ_ONLY, size, nullptr, &operation_result);
    Check(operation_result);
    operation_result = clEnqueueWriteBuffer(queue, buff, CL_FALSE, 0, size, data, 0, nullptr, nullptr);
    Check(operation_result);
}

 
void BHG::upload_force_calc(int num_leaves, int tree_size, int num_bodies){
    using namespace GPU;
    cl_int operation_result = 0 ;
    int mult = 3;
    if (num_leaves == 0) return;

    num_bodies = num_bodies - this -> lower_int_bound - ((this -> upper_int_bound < 0) ? 0 : (num_bodies - this -> upper_int_bound));

    if (warn_if(num_bodies <= 0, "either upper or lower bound were set too high in the gravity solver, resulting in a negative amount of bodies to load")) 
        return; 

    size_t double2_size_tree = sizeof(REAL) * mult * tree_size;
    size_t int_size_tree = sizeof(int) * tree_size;
    size_t double_size_tree = sizeof(REAL) * tree_size;

    size_t double2_size_bodies = sizeof(REAL) * mult * num_bodies; 
    size_t double_size_bodies = sizeof(REAL) * num_bodies;       
    
    cl_mem com_buff, rights_buff, lefts_buff, nmasses_buff, pos_buff, bmasses_buff;

    load_const_buff(double2_size_tree,   this ->nodes.centers_of_mass.data(), com_buff);       
    load_const_buff(int_size_tree,       this ->nodes.left_nodes.data()     , lefts_buff);     
    load_const_buff(int_size_tree,       this ->nodes.right_nodes.data()    , rights_buff);    
    load_const_buff(double_size_tree,    this ->nodes.masses.data()         , nmasses_buff);   
    
    load_const_buff(double2_size_bodies, this -> _s -> bodies.positions.data(),      pos_buff);       
    load_const_buff(double_size_bodies,  this -> _s -> bodies.masses.data(),         bmasses_buff);   


    // accelerations
    cl_mem accs_buff = clCreateBuffer(context,
            CL_MEM_WRITE_ONLY,
            double2_size_bodies,
            nullptr,
        &operation_result);
    Check(operation_result);

    size_t LW_size = 64;
    size_t GW_size = ((num_bodies + LW_size - 1) / LW_size) * LW_size;

    int root_node = num_leaves;  
    REAL tree_bound_size = this ->bounding_box.magnitude(); 

    const int upper = this -> get_upper_bound(), lower = this -> get_lower_bound();
    
    REAL G = this -> _s -> get_G();
    REAL c = this -> _s -> get_c();

    Check(clSetKernelArg(force_calculator, 0, sizeof(int), &num_bodies));
    Check(clSetKernelArg(force_calculator, 1, sizeof(REAL), &G));
    Check(clSetKernelArg(force_calculator, 2, sizeof(REAL), &c));
    Check(clSetKernelArg(force_calculator, 3, sizeof(REAL), &this ->_s ->theta));
    Check(clSetKernelArg(force_calculator, 4, sizeof(REAL), &tree_bound_size));
    Check(clSetKernelArg(force_calculator, 5, sizeof(int), &root_node));      
    Check(clSetKernelArg(force_calculator, 6, sizeof(int), &tree_size));
    Check(clSetKernelArg(force_calculator, 7, sizeof(cl_mem), &nmasses_buff));
    Check(clSetKernelArg(force_calculator, 8, sizeof(cl_mem), &pos_buff)); 
    Check(clSetKernelArg(force_calculator, 9, sizeof(cl_mem), &com_buff)); 
    Check(clSetKernelArg(force_calculator, 10, sizeof(cl_mem), &lefts_buff));
    Check(clSetKernelArg(force_calculator, 11, sizeof(cl_mem), &rights_buff)); 
    Check(clSetKernelArg(force_calculator, 12, sizeof(cl_mem), &bmasses_buff));
    Check(clSetKernelArg(force_calculator, 13, sizeof(cl_mem), &accs_buff));  
    Check(clSetKernelArg(force_calculator, 14, sizeof(int), &lower));
    Check(clSetKernelArg(force_calculator, 15, sizeof(int), &upper)); 


    operation_result = clEnqueueNDRangeKernel(queue, force_calculator, 1, 0, &GW_size, &LW_size, 0, nullptr, nullptr);
    Check(operation_result);


    operation_result = clEnqueueReadBuffer(queue, accs_buff,
        CL_FALSE, 0, 
        double2_size_bodies,
        this ->_s->bodies.accs.data(),
        0, nullptr, nullptr
    ); 
    Check(operation_result); 

    clFinish(queue);
    clReleaseMemObject(com_buff);
    clReleaseMemObject(rights_buff);
    clReleaseMemObject(lefts_buff); 
    clReleaseMemObject(nmasses_buff); 
    clReleaseMemObject(pos_buff); 
    clReleaseMemObject(bmasses_buff);
    clReleaseMemObject(accs_buff);

}

 
void BHG::load_tree_kernel(int n, int num_leaves, void* mortons){
    using namespace GPU;
    cl_int operation_result;
    if (num_leaves == 0) return;
    
    size_t children_size = (n + num_leaves) * sizeof(int);
    size_t mortons_size = num_leaves * sizeof(cl_uint2); 
 
    cl_mem left_nodes_buff  = clCreateBuffer(context, CL_MEM_READ_WRITE, children_size, nullptr, &operation_result);
    Check(operation_result);
    cl_mem right_nodes_buff = clCreateBuffer(context, CL_MEM_READ_WRITE, children_size, nullptr, &operation_result);
    Check(operation_result);
    cl_mem mortons_buff     = clCreateBuffer(context, CL_MEM_READ_ONLY,  mortons_size,  nullptr, &operation_result);
    Check(operation_result);
    
    Check(clEnqueueWriteBuffer(queue, mortons_buff, CL_FALSE, 0, mortons_size, 
                                  mortons, 0, nullptr, nullptr));
    Check(clEnqueueWriteBuffer(queue, left_nodes_buff,  CL_FALSE, 0, children_size,
                                  this ->nodes.left_nodes.data(), 0, nullptr, nullptr));
    Check(clEnqueueWriteBuffer(queue, right_nodes_buff, CL_FALSE, 0, children_size,
                                  this ->nodes.right_nodes.data(), 0, nullptr, nullptr));


    size_t LW_size = 64;
    size_t GW_size = ((n + LW_size - 1) / LW_size) * LW_size;

    operation_result = clSetKernelArg(tree_builder, 0, sizeof(int), &num_leaves);
    Check(operation_result);
    operation_result = clSetKernelArg(tree_builder, 1, sizeof(cl_mem), &mortons_buff);
    Check(operation_result);
    operation_result = clSetKernelArg(tree_builder, 2, sizeof(cl_mem), &left_nodes_buff);
    Check(operation_result);
    operation_result = clSetKernelArg(tree_builder, 3, sizeof(cl_mem), &right_nodes_buff);
    Check(operation_result);


    operation_result = clEnqueueNDRangeKernel(queue, tree_builder, 1, 0, &GW_size, &LW_size, 0, nullptr, nullptr );
    Check(operation_result);

    cl_event evt[2];
    Check(clEnqueueReadBuffer(queue, right_nodes_buff, CL_FALSE, 0, children_size,
                                 this ->nodes.right_nodes.data(), 0, nullptr, &evt[0]));
    Check(clEnqueueReadBuffer(queue, left_nodes_buff,  CL_FALSE, 0, children_size,
                                 this ->nodes.left_nodes.data(),  0, nullptr, &evt[1]));
    
    Check(clWaitForEvents(2, evt));
    clReleaseEvent(evt[0]); clReleaseEvent(evt[1]);
    
    Check(clFinish(queue));
    clReleaseMemObject(left_nodes_buff);
    clReleaseMemObject(right_nodes_buff);
    clReleaseMemObject(mortons_buff);

}


inline void translate2nodes(Simulation* _s, NodesArray& base_layer, const std::vector<cl_uint2>& mortons){
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

 
void BHG::build_tree(){
    this ->threads.clear();
    this ->nodes.clear();
    this ->mortons.clear();

    size_t N = this ->_s->bodies.positions.size();
    if (N == 0) return;
    
    std::vector<cl_uint2> enhanced_mortons;
    enhanced_mortons.resize(N);
    this ->bounding_box = this -> _s ->get_center() * 2;

    parallel_for(uint(0), uint(N), [this, &enhanced_mortons](uint i){
        uint32_t morton_code = get_morton(this, this -> _s -> bodies.positions[i]);
        enhanced_mortons[i].s[0] = morton_code;
        enhanced_mortons[i].s[1] = i;
    });

    parallel_sort(enhanced_mortons.begin(), enhanced_mortons.end(), [](const cl_uint2& a, const cl_uint2& b ){ return a.s[0] < b.s[0];});


    translate2nodes(this -> _s, this ->  nodes, enhanced_mortons);
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

    this ->compressed_mortons_size = num_leaves;
}


BHG::BHG(Simulation* _s) {
    this -> _s = _s;
}

/*
* @brief puts together the force calculation kernel


void BHG::compose_force_kernel(){
    static std::string* force_kernels[4] = {&GPU::newton_cl, &GPU::cl_pn1, &GPU::cl_pn2, &GPU::cl_pn25};

    int index = static_cast<int>(this -> _t);
    if (critical_if(index > 3, "Cannot find suitable kernel for custom force calculating function")) exit(-1);

    std::string kernel_code = *force_kernels[index];
    
    kernel_code += GPU::barnes_force_basic;

    static std::string kernel_names  = "barnes_force";

    GPU::compile_kernel(&kernel_names, &kernel_code, force_calculator);
}*/


void BHG::compose_force_kernel(){
    static std::string* force_kernels[4] = {&GPU::newton_cl3d, &GPU::cl3d_pn1, &GPU::cl3d_pn2, &GPU::cl3d_pn25};

    int index = static_cast<int>(this -> _t);
    if (critical_if(index > 3, "Cannot find suitable kernel for custom force calculating function")){ 
        exit(-1); 
    }

    
    std::string kernel_code = *force_kernels[index];
    kernel_code += GPU::barnes_force_basic_3d;

    static std::string kernel_name = "barnes_force";

    force_calculator = GPU::compile_kernel(&kernel_name, &kernel_code, this->_s->softening);
}


/*
* @brief loads the simulation
*/

void BHG::load(){
    GPU::init_opencl();
    std::string tree_kernel_name = "build_tree";
    std::string force_kernel_name = "barnes_force"; 

    tree_builder = GPU::compile_kernel(&tree_kernel_name, &GPU::barnes_tree_cl, this->_s->softening);
    
    this -> compose_force_kernel();
}


void BHG::compute_forces(){
    this -> build_tree(); 
    upload_force_calc(this -> compressed_mortons_size, this ->nodes.size(), this -> _s -> bodies.positions.size());
}

}   
}