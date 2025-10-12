#pragma once

#include <vector>
#include <cassert>
#include <thread>

#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#define CL_TARGET_OPENCL_VERSION 300

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif

#include "Aster/simulations/barnes-hut.h"

namespace Aster{   
namespace Barnes{


class BH_hyper : public Barnes_Hut {
    public: 

    BH_hyper(Simulation* _s);
    void load() override;
    void compute_forces() override;
    BH_hyper* always_read_positions();
    BH_hyper* read_positions();

    void upload_force_calc(int num_leaves);

    ~BH_hyper();

    private:
    bool always_read_pos = false;
    cl_mem mortons;
    cl_mem sorted_mortons;
    cl_mem left_nodes_buff, right_nodes_buff;
    cl_mem com_buffer, node_masses; 
    cl_mem parents, counters;
    cl_mem bd_size_gpu;


    cl_kernel tree_builder; 

    cl_kernel update_kernel, update_kernel2;
    cl_kernel morton_writer_kernel;
    cl_kernel bottom_up_kernel;  
    cl_kernel first_pass_kernel;
    cl_kernel second_pass_kernel;
    cl_kernel sort_local, sort_global, sort_start;
    cl_kernel force_calculator; 
    
    size_t actual_msize = 1;
    size_t N = 0;

    void compose_sorter();
    void compose_force_kernel();

    void build_bottom_up(int N);
    void load_buffers();
    void make_mortons();
    void sort_mortons();
    void load_tree_kernel(int N);
    void upload_update_kernel(cl_kernel& k, REAL c, REAL d);


    void pull_down_tree(NodesArray& arr);
    void push_tree();

};


}
}