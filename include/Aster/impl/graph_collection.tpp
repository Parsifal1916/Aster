#include <vector>
#include <fstream> 
#include <functional>
#include <cassert>
#include <iostream>
#include <thread>

#include "Aster/graphs/graph_collection.h"
#include "Aster/simulations/sim_obj.h"
#include "Aster/physics/body.h"

namespace Aster{

namespace Graphs{

template <typename T>
Graph<T>::Graph(Simulation<T>* _s, typename Graph<T>::listener_fptr listener, bool for_each_body)
: _s(_s), listener(listener), for_each_body(for_each_body) {
    assert(_s  && "No simulation object specified! seams to be nullptr");
}

template <typename T>
void Graph<T>::init(){
    assert(_s -> bodies.size() && "No bodies to listen to!");
    assert(!done && "already loaded the graph!");

    done = true;

    if (save){
        file = std::ofstream(name + ".csv");
        file << "time,";
    }

    if (for_each_body) {
        size_t num = _s -> bodies.size();
        if (num >= (int)10e3)
            std::cout << "[ ! ] WARN: specified for_each_body=true, the amount of bodies in the simulation results to be very high (>=10e3). large amounts of data are going to be written to disk";
    
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
    assert(data.size() && "empty graphing data! possible data corruption occurred");
    if (for_each_body){ 
        for (int timestep = 0; timestep < data[0].size(); ++timestep){

            file << _s -> get_time_passed() + timestep - data[0].size() << ","; 
            for (int container = 0; container < data.size(); container++)
                file << data[container][timestep] << ((timestep == data[0].size() - 1) ? "\n" : ",");
        }

        for (auto& c : data)
            c.clear();

        return;
    }

    for (int i = 0; i < data[0].size(); i++){
        file << _s -> get_time_passed() + i - data[0].size() << ","; 
        file << data[0][i] << "\n";
    }
    
    data[0].clear();
}

template <typename T>
void Graph<T>::trigger(){
    if (!done) init();
    assert(data.size() && "The internal graphing data is empty!");

    if (data[0].size() >= buffer_size)
        flush_to_file();

    update_data();
}

template <typename T>
void Graph<T>::update_data(){
    if (for_each_body){
        assert(_s -> bodies.size() == data.size() && "The number of bodies has changed over time, cannot generate the graph properly");

        std::vector<std::thread> ts;
        int step = _s -> bodies.size() / _s -> data.NUM_THREADS; 
        int tot = _s -> bodies.size();

        for (int i = 0; i < _s -> data.NUM_THREADS; ++i)
            ts.emplace_back([this, i, step, tot](){
                int start = step*i;
                int stop = (step*(i+2) > tot) ? tot : step*(i*1);

                for (int j = start; j < stop; j++)
                    this -> data[j].push_back(listener(this, this -> _s, &(this -> _s -> bodies[j])));    
            });

        for (auto& t : ts)
            t.join();

        return;
    }

    data[0].push_back(listener(this, _s, nullptr));
}

}
}