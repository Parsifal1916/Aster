#pragma once 
#include <vector>
#include <string>
#include <fstream> 
#include <functional>

#include "Aster/physics/body.h"

namespace Aster{
template <typename T> struct Simulation;

enum graph_type: int {ONCE, FOR_EACH, BETWEEN};

namespace Graphs{

template <typename T>
struct Graph{
    using listener_fptr = std::function<double(struct Graph<T>*, Simulation<T>*  ,Body<T>*)>;
    using collector_fptr = std::function<double(struct Graph<T>*, Simulation<T>*  ,Body<T>*, Body<T>*)>;

    Simulation<T>* _s = nullptr; 
    
    listener_fptr listener;
    collector_fptr collector;
    std::vector<std::vector<double>> data;
    std::string name = "Graph";

    Graph(Simulation<T>* _s, listener_fptr listener,  graph_type type = ONCE);
    Graph(Simulation<T>* _s, collector_fptr listener,  graph_type type = ONCE);

    void trigger();
    void trigger_on(Body<T>* b1, Body<T>* b2);
    void end_batch();
    graph_type get_type() const;

    private:
    graph_type type = ONCE;
    std::ofstream file;
    int buffer_size = 1000;
    double internal_counter = 0;
    bool save = true;
    bool done = false;

    void init();
    void flush_to_file();
    void update_data();
};

template <typename T>
double hamiltonian_collector(Graph<T>* g, Simulation<T>* _s, Body<T>* b);

template <typename T>
double get_total_energy(Simulation<T>* _s);

template <typename T>
double error_collector(Graph<T>* g, Simulation<T>* _s, Body<T>* b);

template <typename T>
double distance_collector(Graph<T>* g, Simulation<T>* _s, Body<T>* b);
}
}

#include "Aster/impl/graph_collection.tpp"