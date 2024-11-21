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

#include "src/physics/vectors.h"
#include "src/physics/body.h"
#include "src/physics/physics.h"

#include "src/sim-setup/sim_helper.h"
#include "src/sim-setup/presets.h"
#include "src/sim-setup/simulation.h"

#include "src/graphics-api/canvas.h"
#include "src/setting/setting.h"