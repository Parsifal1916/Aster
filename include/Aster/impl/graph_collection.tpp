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

/**
* @brief initializes the graph
*/
template <typename T>
void Graph<T>::init(){
    if (warn_if(!_s -> bodies.size(), "No bodies to listen to!"))
        return;
    if (warn_if(done, "already loaded the graph!")) 
        return;

    // saves the loading
    done = true;

    if (save){
        // initilizes the file
        file = std::ofstream(name + ".csv");
        file << "time,";
    }

    // initilizes for a for_eac
    if (type == FOR_EACH) {
        size_t num = _s -> bodies.size();

        // too much time 
        warn_if(num >= (int)10e3, "specified type is FOR_EACH, the amount of bodies in the simulation results to be very high (>=10e3). large amounts of data are going to be written to disk");
    
        // preallocates the memory
        data.resize(num);

        // writes a column for every body
        for (int i = 0; i < num; ++i){
            // initializes the data matrix and reserves some space 
            data[i].reserve(buffer_size);
            if (save) 
                // writes to file
                file << "Body" << i << ", ";
            }

        file << "\n";  
    } else {
        // initilizes for a between or a once
        data.resize(1); // only needs a row

        // reserves the first row with some space
        data[0].reserve(buffer_size);

        // writes to file
        if (save) 
            file << "data\n";
    }
}

/**
* @brief flushes the contents of data to file
*/
template <typename T>
void Graph<T>::flush_to_file(){
    if (err_if(!data.size(), "empty graphing data! possible data corruption occurred"))
        return;

    // flushes a for_each
    if (type == FOR_EACH){ 
        for (int timestep = 0; timestep < data[0].size(); ++timestep){

            // keeps track of the time stamp when the mesurement was taken
            file << _s -> get_time_passed() + timestep - data[0].size() << ",";
            
            // writes the row to file
            for (int container = 0; container < data.size() - 1; container++)
                file << data[container][timestep] << ",";

            // ends the row with a \n
            file << data.back()[timestep] << "\n";
        }

        // cleans the data
        for (auto& c : data)
            c.clear();

        return;
    }

    // writes to file a between or a once
    for (int i = 0; i < data[0].size(); i++){
        // writes the row
        file << _s -> get_time_passed() + i - data[0].size() +1 << ","; 
        
        // ends it with \n
        file << data[0][i] << "\n";
    }

    // cleans the data
    data[0].clear();
}

/**
* @brief triggers the contained graph listener
*/
template <typename T>
void Graph<T>::trigger(){
    if (critical_if(this -> type == BETWEEN, "cannot trigger on graph with type BETWEEN. HINT: try using trigger_on(Body*, Body*)"))
        exit(-1);

    // if it hasn't yet, it initializes itself    
    if (!done) init();

    if (warn_if(!data.size(), "The internal graphing data is empty!"))
        return;

    // checks for a flush
    if (data[0].size() >= buffer_size)
        flush_to_file();

    // updates the data matrix
    update_data();
}


/**
* @brief triggers a BETWEEN graph between two bodies
* @param b1: first body
* @param b2: second body
*/
template <typename T>
void Graph<T>::trigger_on(Body<T>* b1, Body<T>* b2){
    if (err_if(this -> type != BETWEEN, "invalid call to trigger_on: the graph is not of type BETWEEN"))
        return;

    // initilizes if it hasn't
    if (!done){
        this -> init();
        done = true;
    }

    // uses the comulative counter
    this -> internal_counter += collector(this, _s, b1, b2);
}

/**
* @brief same as update_data but for BETWEEN graphs
*/
template <typename T>
void Graph<T>::end_batch(){
    if (err_if(this -> type != BETWEEN, "invalid call to end_batch: the graph is not of type BETWEEN"))
        return;

    // pushes the new batch to the end of the data 
    this -> data[0].push_back(this -> internal_counter);
    internal_counter = 0; // resets the counter
    
    // 
    if (data[0].size() >= buffer_size)
        flush_to_file();
}


/**
* @brief updates the data for a FOR_EACH or for a BETWEEN
*/
template <typename T>
void Graph<T>::update_data(){
    if (err_if(this -> type == BETWEEN, "cannot call update data on graph with type BETWEEN"))
        return;

    if (type == FOR_EACH){
        if (err_if(_s -> bodies.size() != data.size(),"The number of bodies has changed over time, cannot generate the graph properly"))
            return;
        
        // updates the data with a listener for each body
        for (int i = 0; i < _s -> bodies.size(); ++i)
            this -> data[i].push_back(listener(this, this -> _s, &(this -> _s -> bodies[i])));    

        return;
    }

    // calls listener only once
    data[0].push_back(listener(this, _s, nullptr));
}

/**
* @brief returns the type of the graph
* @returns the type of the graph
*/
template <typename T> 
graph_type Graph<T>::get_type() const{
    return this -> type;
}

}
}