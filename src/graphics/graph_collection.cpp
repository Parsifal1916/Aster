#include <vector>
#include <fstream> 
#include <functional>
#include <cassert>
#include <iostream>

#include "Aster/graphs/graph_collection.h"
#include "Aster/simulations/sim_obj.h"
#include "Aster/simulations/3d_sim_obj.h"
#include "Aster/physics/body.h"

namespace Aster{
namespace Graphs{

//using listener2d_fptr = std::function<double(struct Graph2d&, Simulation*  ,Body*)>;
//using listener3d_fptr = std::function<double(struct Graph3d&, Simulation3d*,Body3d*)>;



Graph2d::Graph2d(Simulation* _s, listener2d_fptr listener, bool for_each_body)
: _s(_s), listener(listener), for_each_body(for_each_body) {
    assert(listeners.size() && "Must have at least one listener running! expected >=1 got 0");
    assert(_s  && "No simulation object specified! seams to be nullptr");
    assert(_s -> bodies.size() && "No bodies to listen to!");

    if (save){
        file = std::ofstream(name);
        file << "time,";
    }

    if (for_each_body) {
        size_t num = _s -> bodies.size();
        if (num > (int)10e3)
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

void Graph2d::flush_to_file(){
    assert(data.size() && "empty graphing data! possible data corruption occurred");

    file << _s -> get_time_passed() << ","; 
    if (for_each_body){
        for (auto& container : data){
            for (int i = 0; i < container.size(); i++)
                file << container[i] << (i == container.size() - 1) ? "\n" : ",";

            container.clear();
        }
        return;
    }

    for (int i = 0; i < data[0].size(); i++)
        file << data[0][i] << (i == data[0].size() - 1) ? "\n" : ",";
}

void Graph2d::trigger(){
    assert(data.size() && "The internal graphing data is empty!");

    if (data[0].size() >= buffer_size)
        flush_to_file();

    update_data();
}

void Graph2d::update_data(){
    if (for_each_body){
        assert(bodies.size() == data.size() && "The number of bodies has changed over time, cannot generate the graph properly");

        for (int i = 0; i < _s -> bodies.size(); i++)
            data[i].push_back(listener(this, _s, &(_s -> bodies[i])));

        return;
    }

    data[0].push_back(listener(this, _s, nullptr));
}



}
}