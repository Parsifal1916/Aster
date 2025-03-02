#include <vector>
#include <fstream> 
#include <functional>
#include <cassert>
#include <iostream>
#include <thread>

#include "Aster/graphs/graph_collection.h"
#include "Aster/simulations/sim_obj.h"
#include "Aster/physics/body.h"
#include "Aster/building-api/logging.h"

namespace Aster{

namespace Graphs{

template <typename T>
Graph<T>::Graph(Simulation<T>* _s, typename Graph<T>::listener_fptr listener, graph_type type)
: _s(_s), listener(listener), type(type) {
    critical_if(!_s, "No simulation object specified! seams to be nullptr");
}

template <typename T>
Graph<T>::Graph(Simulation<T>* _s, typename Graph<T>::collector_fptr listener, graph_type type)
: _s(_s), collector(listener), type(type) {
    critical_if(!_s, "No simulation object specified! seams to be nullptr");
}

template <typename T>
void Graph<T>::init(){
    critical_if(!_s -> bodies.size(), "No bodies to listen to!");
    critical_if(done, "already loaded the graph!");

    done = true;

    if (save){
        file = std::ofstream(name + ".csv");
        file << "time,";
    }

    if (type == FOR_EACH) {
        size_t num = _s -> bodies.size();
        warn_if(num >= (int)10e3, "[ ! ] WARN: specified type is FOR_EACH, the amount of bodies in the simulation results to be very high (>=10e3). large amounts of data are going to be written to disk");
    
        data.resize(num);
        for (int i = 0; i < num; ++i){
            data[i].reserve(buffer_size);
            if (save) 
                file << "Body" << i << ", ";
            }

        file << "\n";  
    } else {
        data.resize(1);
        data[0].reserve(buffer_size);
        if (save) 
            file << "data\n";
    }
}

template <typename T>
void Graph<T>::flush_to_file(){
    critical_if(!data.size(), "empty graphing data! possible data corruption occurred");

    if (type == FOR_EACH){ 
        for (int timestep = 0; timestep < data[0].size(); ++timestep){

            file << _s -> get_time_passed() + timestep - data[0].size() << ","; 
            for (int container = 0; container < data.size() - 1; container++)
                file << data[container][timestep] << ",";

            file << data.back()[timestep] << "\n";
        }

        for (auto& c : data)
            c.clear();

        return;
    }

    for (int i = 0; i < data[0].size(); i++){
        file << _s -> get_time_passed() + i - data[0].size() +1 << ","; 
        file << data[0][i] << "\n";
    }

    data[0].clear();
}

template <typename T>
void Graph<T>::trigger(){
    critical_if(this -> type == BETWEEN, "cannot trigger on graph with type BETWEEN. HINT: try using trigger_on(Body*, Body*)");

    if (!done) init();
    critical_if(!data.size(), "The internal graphing data is empty!");

    if (data[0].size() >= buffer_size)
        flush_to_file();

    update_data();
}

template <typename T>
void Graph<T>::trigger_on(Body<T>* b1, Body<T>* b2){
    critical_if(this -> type != BETWEEN, "invalid call to trigger_on: the graph is not of type BETWEEN");
    if (!done){
        this -> init();
        done = true;
    }

    this -> internal_counter += collector(this, _s, b1, b2);
}

template <typename T>
void Graph<T>::end_batch(){
    critical_if(this -> type != BETWEEN, "invalid call to end_batch: the graph is not of type BETWEEN");

    this -> data[0].push_back(this -> internal_counter);
    internal_counter = 0;
    
    if (data[0].size() >= buffer_size)
        flush_to_file();
}



template <typename T>
void Graph<T>::update_data(){
    critical_if(this -> type == BETWEEN, "cannot call update data on graph with type BETWEEN");

    if (type == FOR_EACH){
        critical_if(_s -> bodies.size() != data.size(),"The number of bodies has changed over time, cannot generate the graph properly");

       for (int i = 0; i < _s -> bodies.size(); ++i)
            this -> data[i].push_back(listener(this, this -> _s, &(this -> _s -> bodies[i])));    

        return;
    }

    data[0].push_back(listener(this, _s, nullptr));
}

template <typename T> 
graph_type Graph<T>::get_type() const{
    return this -> type;
}

}
}