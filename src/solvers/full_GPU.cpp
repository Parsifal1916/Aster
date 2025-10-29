#include <type_traits>
#include <algorithm>
#include <unordered_map>
#include <iterator>
#include <future>
#include <chrono>
#include <cmath>

#include <tbb/parallel_for.h>
using namespace tbb;

#include "Aster/simulations/full_GPU.h"

#include "Aster/physics/coeffs.h"
#include "Aster/kernels/BHH_specific.cl.h"
#include "Aster/building-api/GPU_endpoint.h"

#include "Aster/kernels/newton_update.cl.h"
#include "Aster/kernels/newton_update3d.cl.h"
#include "Aster/kernels/PNs.cl.h"
#include "Aster/kernels/PNs3d.cl.h"
#include "Aster/kernels/radix.cl.h"
#include "Aster/kernels/generalized_intgrt.cl.h"
#include "Aster/kernels/generalized_intgrt3d.cl.h"
#include "Aster/kernels/morton_writer.cl.h"

#include <bitset> 

namespace Aster{
namespace Barnes{ 

inline static constexpr auto LOCAL_SIZE = 128;

inline void pull_down(cl_mem buff, void* addr, size_t size){
    using namespace GPU;
    clFinish(queue);
    Check(clEnqueueReadBuffer(queue, buff, CL_TRUE, 0, size, addr, 0, nullptr, nullptr));
}

inline void push_up(cl_mem buff, void* addr, size_t size){
    using namespace GPU;
    Check(clEnqueueWriteBuffer(queue, buff, CL_TRUE, 0, size, addr, 0, nullptr, nullptr));
}



BH_hyper::BH_hyper(Simulation* _s) {
    // sets up the simulation
    this -> _s = _s;
}



void BH_hyper::sort_mortons(){
    using namespace GPU;
    size_t localSize = LOCAL_SIZE;

    size_t globalStart = ( (actual_msize/2 + localSize - 1) / localSize ) * localSize;

    size_t globalGlobal = ( (actual_msize + localSize - 1) / localSize ) * localSize;

    Check(clSetKernelArg(sort_start, 0, sizeof(cl_mem), &mortons));
    Check(clSetKernelArg(sort_start, 1, sizeof(cl_mem), &sorted_mortons));
    Check(clSetKernelArg(sort_start, 2, sizeof(cl_uint), &actual_msize));
    Check(clEnqueueNDRangeKernel(queue, sort_start, 1, NULL, &globalStart, &localSize, 0, NULL, NULL));

    unsigned int limit = 2 * LOCAL_SIZE; 
    for (unsigned int blocksize = limit; blocksize <= actual_msize; blocksize <<= 1) {
        for (unsigned int stride = blocksize / 2; stride > 0; stride >>= 1) {
            Check(clSetKernelArg(sort_global, 0, sizeof(cl_mem), &sorted_mortons));
            Check(clSetKernelArg(sort_global, 1, sizeof(cl_uint), &actual_msize));
            Check(clSetKernelArg(sort_global, 2, sizeof(cl_uint), &blocksize));
            Check(clSetKernelArg(sort_global, 3, sizeof(cl_uint), &stride));
            Check(clEnqueueNDRangeKernel(queue, sort_global, 1, NULL, &globalGlobal, &localSize, 0, NULL, NULL));
        }
    }
}


BH_hyper::~BH_hyper(){
    clReleaseMemObject(mortons);
    clReleaseMemObject(sorted_mortons);
    clReleaseMemObject(node_masses);
    clReleaseMemObject(counters);
    clReleaseMemObject(parents);
    clReleaseMemObject(left_nodes_buff);
    clReleaseMemObject(right_nodes_buff);
    clReleaseMemObject(com_buffer);   
}

 
void BH_hyper::make_mortons(){
    using namespace GPU;

    int N = this  -> _s->bodies.positions.size();
    static int prec_bits = PRECISION_BITS;

    const int lower = this->get_lower_bound(), upper = this->get_upper_bound();

    size_t localSize  = LOCAL_SIZE;
    size_t globalSize = ((N + localSize - 1) / localSize) * localSize;
    Check(clSetKernelArg(morton_writer_kernel, 0, sizeof(uint), &N));
    Check(clSetKernelArg(morton_writer_kernel, 1, sizeof(cl_mem), &this -> _s -> positions_cl));
    Check(clSetKernelArg(morton_writer_kernel, 2, sizeof(cl_mem), &mortons));
    Check(clSetKernelArg(morton_writer_kernel, 3, sizeof(cl_mem), &bd_size_gpu));
    Check(clSetKernelArg(morton_writer_kernel, 4, sizeof(uint), &lower));
    Check(clSetKernelArg(morton_writer_kernel, 5, sizeof(uint), &upper));
    
    Check(clEnqueueNDRangeKernel(queue, morton_writer_kernel, 1, NULL, &globalSize, &localSize, 0, NULL, NULL));
    clFinish(queue);
}

/*
* @brief puts together the force calculation kernel

void BH_hyper<vec2>::compose_force_kernel(){
    static std::string* force_kernels[4] = {&GPU::newton_cl, &GPU::cl_pn1, &GPU::cl_pn2, &GPU::cl_pn25};

    int index = static_cast<int>(this -> _t);
    if (critical_if(index > 3, "Cannot find suitable kernel for custom force calculating function")) exit(-1);

    std::string kernel_code = *force_kernels[index];
    
    kernel_code += GPU::BHH_force_basic;

    static std::string kernel_names  = "barnes_force";

    GPU::compile_kernel(&kernel_names, &kernel_code, force_calculator);
}
*/


void BH_hyper::compose_force_kernel(){
    static std::string* force_kernels[4] = {&GPU::newton_cl3d, &GPU::cl3d_pn1, &GPU::cl3d_pn2, &GPU::cl3d_pn25};

    int index = static_cast<int>(this -> _t);
    if (critical_if(index > 3, "Cannot find suitable kernel for custom force calculating function")) exit(-1);

    
    std::string kernel_code = *force_kernels[index];
    kernel_code += GPU::BHH_force_basic3d;

    static std::string kernel_name = "barnes_force";

    force_calculator = GPU::compile_kernel(&kernel_name, &kernel_code, this->_s->softening);
}


void BH_hyper::compose_sorter(){
    using namespace GPU;

    static std::string m = "gen_mortons";
    morton_writer_kernel = compile_kernel(&m, &morton_writer3d, this->_s->softening);

    static std::string fpk = "first_pass";
    first_pass_kernel = compile_kernel(&fpk, &bottom_up3d, this->_s->softening);
 
    static std::string bs = "Sort_BitonicMergesortStart";
    sort_start = compile_kernel(&bs, &radix_cl, this->_s->softening);

    static std::string bg = "Sort_BitonicMergesortGlobal";
    sort_global = compile_kernel(&bg, &radix_cl, this->_s->softening);

    std::string tree_kernel_name = "build_tree";
    tree_builder = compile_kernel(&tree_kernel_name, &GPU::BHH_tree_cl3d, this->_s->softening);
}


void BH_hyper::build_bottom_up(int N){
    using namespace GPU;
    size_t LW_size = 64;
    size_t leaf_only_GW = ((N + LW_size - 1) / LW_size) * LW_size;
    
    size_t num = N;
    Check(clSetKernelArg(first_pass_kernel, 0, sizeof(int),    &num));  
    Check(clSetKernelArg(first_pass_kernel, 1, sizeof(cl_mem), &left_nodes_buff));
    Check(clSetKernelArg(first_pass_kernel, 2, sizeof(cl_mem), &right_nodes_buff));
    Check(clSetKernelArg(first_pass_kernel, 3, sizeof(cl_mem), &com_buffer));
    Check(clSetKernelArg(first_pass_kernel, 4, sizeof(cl_mem), &node_masses));
    Check(clSetKernelArg(first_pass_kernel, 5, sizeof(cl_mem), &counters));
    Check(clSetKernelArg(first_pass_kernel, 6, sizeof(cl_mem), &parents));


    Check(clEnqueueNDRangeKernel(queue, first_pass_kernel, 1, 0, &leaf_only_GW, &LW_size, 0, nullptr, nullptr));
    clFinish(queue);
}

void compose_updater(){
    using namespace GPU;
    static std::string update_kernel_src = general_saba;

    static std::string kernel_name = "saba";
    int index =  static_cast<int>(1);


    if (critical_if(index < 0 || index >= static_cast<int>(CUSTOM_U), 
    "could not find a suitable GPU-accelerated function for kernel " + kernel_name))
        exit(1);

    //compile_kernel(&kernel_name, &update_kernel_src, update_kernel);
}
/*
void compose_updater3(){
    using namespace GPU;
    static std::string update_kernel_src = general_saba3d;

    static std::string kernel_name = "saba";
    int index =  static_cast<int>(1);


    if (critical_if(index < 0 || index >= static_cast<int>(CUSTOM_U), 
    "could not find a suitable GPU-accelerated function for kernel " + kernel_name))
        exit(1);

    //compile_kernel(&kernel_name, &update_kernel_src, update_kernel);
}*/


void BH_hyper::load_buffers(){
    using namespace GPU;
    size_t N = this  -> _s -> bodies.positions.size();
    int mult = 3;      

    int vec_array = N * mult * sizeof(REAL);
    int mrt_array = actual_msize * sizeof(cl_ulong);
    int child_array = (2*N-1) * sizeof(cl_int);
    int double_array = sizeof(REAL) * N;


    cl_int err;
    
    //----------------------------  MORTON SORTING ------------------------------
    std::vector<uint64_t> maxes(actual_msize, 0 );
    sorted_mortons  = clCreateBuffer(context, CL_MEM_READ_WRITE, mrt_array, NULL, &err);
    Check(err);
    mortons         = clCreateBuffer(context, CL_MEM_READ_WRITE, mrt_array, NULL, &err);
    Check(err);
    parents         = clCreateBuffer(context, CL_MEM_READ_WRITE, (2*N -1) * sizeof(cl_uint), NULL, &err);
    Check(err);
    bd_size_gpu     = clCreateBuffer(context, CL_MEM_READ_WRITE, 4 * sizeof(REAL), NULL, &err);
    Check(err);
    Check(clEnqueueWriteBuffer(queue, bd_size_gpu,  CL_TRUE, 0, sizeof(vec3), &this -> bounding_box, 0, NULL, NULL));

    counters         =  clCreateBuffer(context, CL_MEM_READ_WRITE, (2*N -1) * sizeof(cl_uint), NULL, &err);
    Check(err);
    

    //---------------------------- LEFT/RIGHT NODE PTRS ----------------------------
    left_nodes_buff  = clCreateBuffer(context, CL_MEM_READ_WRITE, child_array, nullptr, &err);
    Check(err);
    right_nodes_buff = clCreateBuffer(context, CL_MEM_READ_WRITE, child_array, nullptr, &err);
    Check(err);

    std::vector<int> neg_one(2*N-1, -1);
    clEnqueueWriteBuffer(queue, right_nodes_buff, CL_TRUE, 0, child_array, neg_one.data(), 0, NULL, NULL);
    clEnqueueWriteBuffer(queue, left_nodes_buff, CL_TRUE, 0,  child_array, neg_one.data(), 0, NULL, NULL);

    //--------------------------- MORTON TREE SOA -----------------------------------
    com_buffer      = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(vec3) * (2*N -1), NULL, &err);
    Check(err);
    node_masses     = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(REAL) * (2*N -1), NULL, &err);
    Check(err);
}
/*
* @brief loads the simulation
*/

void BH_hyper::load(){
    this -> bounding_box = this -> _s -> get_center()*2;

    // assigns base value for the total mass
    this -> _s -> calculate_total_mass();

    GPU::init_opencl();
    this ->_s-> load_gpu_buffers();
    
    this -> N = this  -> _s -> bodies.positions.size();
    this -> nodes.resize(N*2 - 1);
    
    actual_msize = 128;
    while (actual_msize < N) actual_msize *=2;

    this -> compose_sorter();
    //this -> compose_updater();
    this -> compose_force_kernel();

    this -> load_buffers();
}

void BH_hyper::upload_force_calc(int num_leaves){
    using namespace GPU;
    cl_int operation_result = 0 ;
    int num_bodies = num_leaves;
    REAL dt = this -> _s -> get_dt();

    size_t LW_size = 64;
    size_t GW_size = ((N + LW_size - 1) / LW_size) * LW_size;

    REAL init_size = this->bounding_box.magnitude();  
    REAL theta = this->_s->theta;
    REAL G = this -> _s -> get_G();
    REAL c = this -> _s -> get_c();

    int lower = this->get_lower_bound(), upper = this->get_upper_bound();
    
    Check(clSetKernelArg(force_calculator, 0,  sizeof(int),    &num_bodies       ));
    Check(clSetKernelArg(force_calculator, 1,  sizeof(REAL),   &init_size        ));  
    Check(clSetKernelArg(force_calculator, 2,  sizeof(REAL),   &G                ));
    Check(clSetKernelArg(force_calculator, 3,  sizeof(REAL),   &c                ));
    Check(clSetKernelArg(force_calculator, 4,  sizeof(REAL),   &theta            ));
    Check(clSetKernelArg(force_calculator, 5,  sizeof(cl_mem), &left_nodes_buff  ));
    Check(clSetKernelArg(force_calculator, 6,  sizeof(cl_mem), &right_nodes_buff )); 
    Check(clSetKernelArg(force_calculator, 7,  sizeof(cl_mem), &node_masses      ));
    Check(clSetKernelArg(force_calculator, 8,  sizeof(cl_mem), &com_buffer       )); 
    Check(clSetKernelArg(force_calculator, 9,  sizeof(cl_mem), &this -> _s -> positions_cl)); 
    Check(clSetKernelArg(force_calculator, 10, sizeof(cl_mem), &this -> _s -> masses_cl ));
    Check(clSetKernelArg(force_calculator, 11, sizeof(cl_mem), &this -> _s -> accs_cl  ));  
    Check(clSetKernelArg(force_calculator, 12,  sizeof(int),   &lower            ));
    Check(clSetKernelArg(force_calculator, 13,  sizeof(int),   &upper            ));

    Check(clEnqueueNDRangeKernel(queue, force_calculator, 1, 0, &GW_size, &LW_size, 0, nullptr, nullptr));
    clFinish(queue);
}

 
void BH_hyper::load_tree_kernel(int N){
    using namespace GPU;
    cl_int operation_result;
    int num_leaves = N;
    if (!num_leaves) return;


    size_t LW_size = 64;
    uint num_internal = num_leaves - 1;
    size_t GW_size = ((N + LW_size - 1) / LW_size) * LW_size;

    Check(clSetKernelArg(tree_builder, 0, sizeof(cl_int), &num_leaves));
    Check(clSetKernelArg(tree_builder, 1, sizeof(cl_mem), &sorted_mortons));
    Check(clSetKernelArg(tree_builder, 2, sizeof(cl_mem), &left_nodes_buff));
    Check(clSetKernelArg(tree_builder, 3, sizeof(cl_mem), &right_nodes_buff));
    Check(clSetKernelArg(tree_builder, 4, sizeof(cl_mem),  &com_buffer));
    Check(clSetKernelArg(tree_builder, 5, sizeof(cl_mem),  &node_masses));
    Check(clSetKernelArg(tree_builder, 6, sizeof(cl_mem),  &this -> _s -> positions_cl));
    Check(clSetKernelArg(tree_builder, 7, sizeof(cl_mem),  &this -> _s -> masses_cl));
    Check(clSetKernelArg(tree_builder, 8, sizeof(cl_mem),  &parents));
    Check(clSetKernelArg(tree_builder, 9, sizeof(cl_mem),  &counters));

    operation_result = clEnqueueNDRangeKernel(queue, tree_builder, 1, 0, &GW_size, &LW_size, 0, nullptr, nullptr );
    Check(operation_result);
    clFinish(queue);
}


BH_hyper* BH_hyper::read_positions(){
    using namespace GPU;
    clFinish(queue);
    Check(clEnqueueReadBuffer(queue, this -> _s -> positions_cl, CL_FALSE, 0, sizeof(vec3) * N, this  -> _s -> bodies.positions.data(), 0, nullptr, nullptr));
    return this;
}


BH_hyper* BH_hyper::always_read_positions(){
    always_read_pos = true;
    return this;
}

/*
void BH_hyper::debug_step() {
    using namespace GPU;
    using clock = std::chrono::high_resolution_clock;
    using ms = std::chrono::duration<double, std::milli>;


    static const REAL* sequence = saba_coeffs[static_cast<int>(this->update_used)];
    int N = this->bodies.positions.size();
    int order = static_cast<int>(this->update_used);
    static int mult = sizeof(T);

    std::vector<std::pair<std::string, double>> timings;

    auto time_block = [&](const std::string& name, auto&& func) {
        auto start = clock::now();
        func();
        clFinish(queue);
        auto end = clock::now();
        timings.emplace_back(name, ms(end - start).count());
    };

    time_block("make_mortons", [&]() { make_mortons();});
    time_block("sort_mortons", [&]() { sort_mortons(); });
    time_block("load_tree_kernel", [&]() { load_tree_kernel(N); }); 
    time_block("upload_force_calc", [&]() {upload_force_calc(N); });
    time_block("update_bodies", [&]() {upload_update_kernel(update_kernel, 1, 1);});

    if (always_read_pos) {
        time_block("read_positions", [&]() {
            Check(clEnqueueReadBuffer(queue, positions, CL_FALSE, 0,
                                      sizeof(T) * N, this->bodies.positions.data(),
                                      0, nullptr, nullptr));
        });
    }

    std::sort(timings.begin(), timings.end(),
              [](auto& a, auto& b) { return a.second > b.second; });

    std::cout << "\n=== Profiling step report ===\n";
    for (auto& [name, time] : timings) {
        std::cout << name << ": " << time << " ms\n";
    }
    std::cout << "=============================\n";
}
*/


void BH_hyper::pull_down_tree(NodesArray& arr){
    this ->nodes.resize(2*N-1);
    pull_down(left_nodes_buff,  arr.left_nodes.data(),      (2*N-1) * sizeof(int));
    pull_down(right_nodes_buff, arr.right_nodes.data(),     (2*N-1) * sizeof(int));
    pull_down(node_masses,      arr.masses.data(),          (2*N-1) * sizeof(REAL));
    pull_down(com_buffer,       arr.centers_of_mass.data(), (2*N-1) * sizeof(vec2));
}


void BH_hyper::push_tree(){
    this ->nodes.resize(2*N-1);
    push_up(left_nodes_buff,    this -> nodes.left_nodes.data(),     (2*N-1) * sizeof(cl_int));
    push_up(right_nodes_buff,   this -> nodes.right_nodes.data(),    (2*N-1) * sizeof(cl_int));
    push_up(node_masses,        this -> nodes.masses.data(),         (2*N-1) * sizeof(REAL));
    push_up(com_buffer,         this -> nodes.centers_of_mass.data(),(2*N-1) * sizeof(vec2));
}


inline void print_tree(NodesArray& tree,NodesArray& tree2, std::string label){
    std::cout << "+++++++++++++ start of " << label << " +++++++++++++\n";
    for (int i = (tree.size() + 1) / 2; i < tree.size(); ++i){
        std::cout << tree.masses[i] << "kg " << tree.centers_of_mass[i] << " l" << tree.left_nodes[i] << "r" << tree.right_nodes[i] << "\n";
        //if (tree2.masses[i]  - tree.masses[i] != 0) std::cout << "mismatch " << tree2.masses[i]  - tree.masses[i]<< "\n";
        //if (tree2.centers_of_mass[i].magnitude()  - tree.centers_of_mass[i].magnitude() != 0) std::cout << "mismatch " << tree2.centers_of_mass[i]  - tree.centers_of_mass[i] << "\n";
    }
    std::cout << "+++++++++++++ end of " << label << " +++++++++++++++\n";
}

//
//uint32_t get_tbreaked_morton(Barnes_Hut* _s, const T& pos, size_t body_index) {
//    return get_morton(_s, pos);
//}


inline void t2nodes(Simulation* _s, NodesArray& base_layer, const std::vector<std::pair<uint32_t, unsigned int>>& mortons){
    base_layer.clear();
    base_layer.resize(mortons.size());

    base_layer.init(_s,mortons[0].first, 0);
    size_t last = 0;

    for (size_t i = 1; i < mortons.size(); ++i) {
        last++;
        base_layer.init(_s, mortons[i].first, last);
    }
    base_layer.resize(last+1);
}


inline void cpu_tree(std::vector<std::pair<uint32_t, unsigned int>>& enhanced_mortons, BH_hyper* _s, int N, NodesArray& nodes){
    size_t num_leaves = N;
    size_t num_internal = num_leaves - 1;
    nodes.resize(num_leaves + num_internal);

    auto delta = [&](size_t i, size_t j) -> int {
        if (j >= enhanced_mortons.size()) return -1;
        if (enhanced_mortons[i].first == enhanced_mortons[j].first) {
            return 32 + __builtin_clz(i ^ j); 
        }
        return __builtin_clz(enhanced_mortons[i].first ^ enhanced_mortons[j].first);
    };

    parallel_for(size_t(0), num_internal, [&, delta, _s](size_t i){
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
            _s -> nodes.left_nodes[internal_node] = gamma; 
        } else {
            _s -> nodes.left_nodes[internal_node] = num_leaves + gamma; 
        }
        
        if (right_range == gamma + 1) {
            _s -> nodes.right_nodes[internal_node] = gamma + 1;  
        } else {
            _s -> nodes.right_nodes[internal_node] = num_leaves + gamma + 1; 
        }
    });
}


void BH_hyper::compute_forces() {
    using namespace GPU;
    make_mortons();
    sort_mortons();

    load_tree_kernel(this -> get_range());
    build_bottom_up(this -> get_range());
    upload_force_calc(this -> _s->bodies.positions.size());
}
}   
}

/*
for (size_t i = 0; i < N; ++i) {
        uint32_t morton_code = get_tbreaked_morton(this, this->bodies.positions[i], i);
        enhanced_mortons.emplace_back(i, morton_code);
    }

    std::sort(enhanced_mortons.begin(), enhanced_mortons.end(), [](const std::pair<uint32_t, size_t>& a, const std::pair<uint32_t, size_t>&b ){ return a.second < b.second;});


    t2nodes(this, this -> nodes, enhanced_mortons);
    
    if (N <= 1) return;
  
    size_t num_internal = N - 1;
    this ->nodes.resize(2* N-1);

    auto delta = [&](size_t i, size_t j) -> int {
        if (j >= enhanced_mortons.size()) return -1;
        if (enhanced_mortons[i].second == enhanced_mortons[j].second) {
            return 32 + __builtin_clz(i ^ j); 
        }
        return __builtin_clz(enhanced_mortons[i].second ^ enhanced_mortons[j].second);
    };

    parallel_for(size_t(0), num_internal, [&, delta, this](size_t i){
        int d_left = (i == 0) ? -1 : delta(i, i - 1);
        int d_right = delta(i, i + 1);
        int d = (d_right > d_left) ? 1 : -1;
        
        int d_min = (d == 1) ? d_left : d_right;
        size_t l_max = 2;
        while ((int)(i + l_max * d) >= 0 && i + l_max * d < N && delta(i, i + l_max * d) > d_min) {
            l_max *= 2;
        }
        
        size_t l = 0;
        for (size_t t = l_max / 2; t >= 1; t /= 2) {
            size_t test_pos = i + (l + t) * d;
            if ((int)test_pos >= 0 && test_pos < N && delta(i, test_pos) > d_min) {
                l += t;
            }
        }
        
        size_t j = i + l * d;
        int d_node = delta(i, j);
        
        size_t s = 0;
        size_t range_size = (j > i) ? j - i : i - j;
        for (size_t t = (range_size + 1) / 2; t >= 1; t = (t + 1) / 2) {
            size_t test_pos = i + (s + t) * d;
            if ((int)test_pos >= 0 && test_pos < N && delta(i, test_pos) > d_node) {
                s += t;
            }
            if (t == 1) break;
        }
        
        size_t gamma = i + s * d + std::min(d, 0);
        
        size_t internal_node = N + i;
        
        size_t left_range = std::min(i, j);
        size_t right_range = std::max(i, j);
         
        if (left_range == gamma) {
            this -> nodes.left_nodes[internal_node] = gamma; 
        } else {
            this -> nodes.left_nodes[internal_node] = N + gamma; 
        }
        
        if (right_range == gamma + 1) {
            this ->nodes.right_nodes[internal_node] = gamma + 1;  
        } else {
            this -> nodes.right_nodes[internal_node] = N + gamma + 1; 
        }
    });

    this -> nodes.unite(N);
    this->compressed_mortons_size = N;
*/