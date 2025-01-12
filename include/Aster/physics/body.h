#pragma once

#include "Aster/physics/vectors.h"

namespace Aster{

struct Body{
public:
       double mass = 1;
       double temp = 10;
       vec2 position = {};
       vec2 velocity = {};
       vec2 acceleration = {};
       vec2 prev_acc = {};

       Body(double mass, vec2 position, vec2 velocity) : mass(mass),position(position), velocity(velocity){};
       Body(){}
};


struct Body3d{
       double mass = 1;
       double temp = 10;
       vec3 position = {};
       vec3 velocity = {};
       vec3 acceleration = {};
       vec3 prev_acc = {};

       Body3d(double mass, vec3 position, vec3 velocity) : mass(mass),position(position), velocity(velocity){};
       Body3d(){}

       static float* get_color(double t) {
              static float color[3];

              if (t <= 0.3) {
                     double normalizedT = t / 0.33;
                     color[0] = 0 + normalizedT * (255 - 0);
                     color[1] = 0;
                     color[2] = 255 - normalizedT * (255 - 0);
              } else{
                     double normalizedT = (t - 0.33) / 0.33;
                     color[0] = 255;
                     color[1] = 0 + normalizedT * (255 - 0);
                     color[2] = 0;
              }
              return color;
       }

       static double sigmoid(double x) {
              return 1.0 / (1.0 + exp(-std::log10(x)));
       }
};

}