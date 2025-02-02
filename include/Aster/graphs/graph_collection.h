#pragma once 
#include <vector>
#include <string>
#include <fstream> 
#include <functional>

#include "Aster/physics/body.h"
<<<<<<< HEAD

namespace Aster{
template <typename T> struct Simulation;

namespace Graphs{

template <typename T>
struct Graph{
    using listener_fptr = std::function<double(struct Graph<T>*, Simulation<T>*  ,Body<T>*)>;
    bool for_each_body = false;

    Simulation<T>* _s = nullptr; 
    
    listener_fptr listener;
    std::vector<std::vector<double>> data;
    std::string name = "Graph";

    Graph(Simulation<T>* _s, listener_fptr listener,  bool for_each_body = false);
=======
#include "Aster/simulations/sim_obj.h"
#include "Aster/simulations/3d_sim_obj.h"

namespace Aster{
namespace Graphs{


using listener3d_fptr = std::function<double(struct Graph3d*, Simulation3d*,Body3d*)>;

struct Graph2d{
    using listener2d_fptr = std::function<double(struct Graph2d*, Simulation*  ,Body*)>;
    bool for_each_body = false;

    Simulation* _s = nullptr; 
    
    listener2d_fptr listener;
    std::vector<std::vector<double>> data;
    std::string name = "Graph";

    Graph2d(Simulation* _s, listener2d_fptr listener,  bool for_each_body = false);
>>>>>>> main

    void trigger();

    private:
    std::ofstream file;
    int buffer_size = 1000;
    bool save = true;
<<<<<<< HEAD
    bool done = false;

    void init();
=======

>>>>>>> main
    void flush_to_file();
    void update_data();
};

}
<<<<<<< HEAD
}

#include "Aster/impl/graph_collection.tpp"
=======
}
>>>>>>> main
