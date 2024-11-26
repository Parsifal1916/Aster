// *  █████╗ ███████╗████████╗███████╗██████╗ 
// * ██╔══██╗██╔════╝╚══██╔══╝██╔════╝██╔══██╗
// * ███████║███████╗   ██║   █████╗  ██████╔╝
// * ██╔══██║╚════██║   ██║   ██╔══╝  ██╔══██╗
// * ██║  ██║███████║   ██║   ███████╗██║  ██║
// * ╚═╝  ╚═╝╚══════╝   ╚═╝   ╚══════╝╚═╝  ╚═╝

#include <SFML/Graphics.hpp>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <chrono>
#include <iomanip>
#include <thread>
#include <ctime>

#include "physics/vectors.h"
#include "physics/body.h"
#include "physics/physics.h"

#include "sim-setup/sim_helper.h"
#include "sim-setup/presets.h"
#include "sim-setup/simulation.h"

#include "graphics-api/canvas.h"
#include "setting/setting.h"

int main(){
    Simulation::init();
    Simulation::presets::add_body(5e11, {WIDTH/8 + WIDTH/9, HEIGHT/2}, {0,  200}, sf::Color::Red);
    Simulation::presets::add_body(5e11, {WIDTH/7 + WIDTH/9, HEIGHT/2}, {0, -200}, sf::Color::Blue);
    Simulation::fire_up();
    
}