#pragma once

#include <cmath>
#include <vector>
#include "Aster/physics/vectors.h"

namespace Aster{

struct BodyArray{
public:
       std::vector<REAL> masses = {};
       std::vector<REAL> temps = {};
       std::vector<vec3> positions = {};
       std::vector<vec3> velocities = {};
       std::vector<vec3> accs = {};

       BodyArray(){}

       vec3& get_position_of(size_t i);
       vec3& get_velocity_of(size_t i);
       REAL& get_mass_of(size_t i);
       REAL& get_temp_of(size_t i);
       vec3& get_acc_of(size_t i);
};

}