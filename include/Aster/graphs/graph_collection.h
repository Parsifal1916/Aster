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
    // func pointer to a ONCE or FOR_EACH  
    using listener_fptr = std::function<double(struct Graph<T>*, Simulation<T>*  ,size_t)>;
    
    // func pointer to a BETWEEN
    using collector_fptr = std::function<double(struct Graph<T>*, Simulation<T>*  ,size_t, size_t)>;

    Simulation<T>* _s = nullptr; 
    
    listener_fptr listener;
    collector_fptr collector;
    std::vector<std::vector<double>> data;
    std::string name = "Graph";

    Graph(Simulation<T>* _s, listener_fptr listener,  graph_type type = ONCE);
    Graph(Simulation<T>* _s, collector_fptr listener,  graph_type type = ONCE);

    /**
    * @brief triggers the contained graph listener
    */
    void trigger();

    /**
    * @brief triggers a BETWEEN graph between two bodies
    * @param b1: first body
    * @param b2: second body
    */
    void trigger_on(size_t b1, size_t b2);

    /**
    * @brief same as update_data but for BETWEEN graphs
    */
    void end_batch();
    graph_type get_type() const;

    private:
    // default type
    graph_type type = ONCE;
    // file to write to
    std::ofstream file;
    // maximum buffer size
    int buffer_size = 50;
    // internal cumulative counter
    double internal_counter = 0;

    bool save = true; // should it save?
    bool done = false; // has it loaded?

    /**
    * @brief initializes the graph
    */
    void init();

    /**
    * @brief flushes the contents of data to file
    */
    void flush_to_file();

    /**
    * @brief updates the data for a FOR_EACH or for a BETWEEN
    */
    void update_data();
};

template <typename T>
double hamiltonian_collector(Graph<T>* g, Simulation<T>* _s, size_t b);

template <typename T>
double get_total_energy(Simulation<T>* _s);

template <typename T>
double error_collector(Graph<T>* g, Simulation<T>* _s, size_t b);

template <typename T>
double distance_collector(Graph<T>* g, Simulation<T>* _s, size_t b);
}
}

#include "Aster/impl/graph_collection.tpp"