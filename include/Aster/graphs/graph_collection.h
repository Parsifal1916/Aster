#pragma once 
#include <vector>
#include <string>
#include <fstream> 
#include <functional>

#include "Aster/physics/body.h"

namespace Aster{
struct Simulation;
struct Simulation3d;

namespace Graphs{

struct Graph2d{
    using listener2d_fptr = std::function<double(struct Graph2d*, Simulation*  ,Body*)>;
    bool for_each_body = false;

    Simulation* _s = nullptr; 
    
    listener2d_fptr listener;
    std::vector<std::vector<double>> data;
    std::string name = "Graph";

    Graph2d(Simulation* _s, listener2d_fptr listener,  bool for_each_body = false);

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

struct Graph3d{
    using listener3d_fptr = std::function<double(struct Graph3d*, Simulation3d*  ,Body3d*)>;
    bool for_each_body = false;

    Simulation3d* _s = nullptr; 
    
    listener3d_fptr listener;
    std::vector<std::vector<double>> data;
    std::string name = "Graph";

    Graph3d(Simulation3d* _s, listener3d_fptr listener,  bool for_each_body = false);

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