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

double time_passed = 0;

#define spin .9// v -> 1/spin
#define stV 0

#include "physics/vectors.h"
#include "physics/body.h"
#include "physics/physics.h"

#include "sim-setup/sim_helper.h"
#include "sim-setup/presets.h"
#include "sim-setup/simulation.h"

#include "graphics-api/canvas.h"
#include "setting/setting.h"


int main() {
    Simulation::init();
    Simulation::populate();
    Simulation::fire_up();

    return 0;
}
