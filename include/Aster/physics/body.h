#pragma once

#include <cmath>
#include <vector>
#include "Aster/physics/vectors.h"

namespace Aster{

template <typename T>
struct BodyArray{
public:
       std::vector<double> masses = {};
       std::vector<double> temps = {};
       std::vector<T> positions = {};
       std::vector<T> velocities = {};
       std::vector<T> accs = {};

       BodyArray(){}

       T& get_position_of(size_t i);
       T& get_velocity_of(size_t i);
       double& get_mass_of(size_t i);
       double& get_temp_of(size_t i);
       T& get_acc_of(size_t i);
};

template<typename T>
T& BodyArray<T>::get_position_of(size_t i){
       return positions[i];
}
template<typename T>
T& BodyArray<T>::get_velocity_of(size_t i){
       return velocities[i];
}

template<typename T>
double& BodyArray<T>::get_mass_of(size_t i){
       return masses[i];
}

template<typename T>
double& BodyArray<T>::get_temp_of(size_t i){
       return temps[i];
}

template<typename T>
T& BodyArray<T>::get_acc_of(size_t i){
       return accs[i];
}



}