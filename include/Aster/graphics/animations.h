#pragma once

#include <GLFW/glfw3.h>
#include <thread>
#include <iostream>

#include "Aster/simulations/sim_obj.h"

namespace Aster{
namespace Renderer{

void draw_quad(GLFWwindow* window, int loaded_bodies, int total_bodies);
void show_loadingbar(GLFWwindow* window, Simulation<vec2>* _s );

}
}