#pragma once

#include <cmath>
#include "Aster/physics/vectors.h"

namespace Aster{

template <typename T>
struct Body{
public:
       double mass = 1;
       double temp = 10;
       T position = {};
       T velocity = {};
       T acceleration = {};
       T prev_acc = {};

       Body(double mass, T position, T velocity, double temp = 0) : mass(mass),position(position), velocity(velocity){};
       Body(){}
};

}