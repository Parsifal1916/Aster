#pragma once

#include <GLFW/glfw3.h>
#include <thread>
#include <iostream>

#include "Aster/simulations/sim_obj.h"

namespace Aster{
namespace Renderer{

/**
* @brief draws a loading bar onto a window with a certain percentage
* @param window: window pointer to draw in
* @param percent: loading bar percentage
*/
void draw_quad(GLFWwindow* window, int loaded_bodies, int total_bodies);

/**
* @brief shows a loading bar while loading a simulation's queue
* @param window: window ptr to load in
* @param _s: simulation to load
*/
void show_loadingbar(GLFWwindow* window, Simulation* _s );

}
}