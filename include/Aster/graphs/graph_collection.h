#pragma once 
#include <vector>
#include <string>
#include <fstream> 
#include <functional>

#include "Aster/physics/body.h"

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

    void trigger();

    private:
    std::ofstream file;
    int buffer_size = 1000;
    bool save = true;
    bool done = false;

    void init();
    void flush_to_file();
    void update_data();
};

}
}

#include "Aster/impl/graph_collection.tpp"