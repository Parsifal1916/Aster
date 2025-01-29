#pragma once 
#include <cmath>
#include <vector>
#include <random>
#include <iomanip>
#include <sstream>

#include "Aster/simulations/sim_obj.h"
#include "Aster/physics/body.h"

namespace Aster{

//===---------------------------------------------------------===//
// Backend initialization                                        //
//===---------------------------------------------------------===//

bool 
    leave_traces = false,
    save = false
;

using func_ptr = void(*)(Body*, Simulation*);
using force_func = vec2(*)(double, double, vec2, vec2, vec2, vec2, Simulation*);

std::mt19937 rng;

int 
    WIDTH = 1366, 
    HEIGHT = 768
;

std::string format_time(int totalSeconds) {
    int hours = totalSeconds / 3600;
    int minutes = (totalSeconds % 3600) / 60;
    int seconds = totalSeconds % 60;

    std::ostringstream oss;
    oss << std::setw(2) << std::setfill('0') << hours << ":"
        << std::setw(2) << std::setfill('0') << minutes << ":"
        << std::setw(2) << std::setfill('0') << seconds;
    return oss.str();
}
}