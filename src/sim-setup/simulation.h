#pragma once 
#include <SFML/Graphics.hpp>
#include <random>
#include <cmath>
#include <vector>
#include "../physics/body.h"

struct Body;

namespace Simulation{
    double lagrangian = 0;
    double highest_lagrangian = 0;
    std::vector<double> lagrangians = {0};

    bool leave_traces = false;
    bool save = false;

    auto bg = sf::Color(0,0,0);
    std::vector<std::thread> threads;
    std::vector<Body> bodies;

    enum simulation_types: int {LIGHT = 0, MEDIUM = 1, HEAVY = 2}; 
    simulation_types sim_type = simulation_types::LIGHT;

    using step_func = void(*)();
    using func_ptr = void(*)(Body*);
    using force_func = vec2(*)(double, double, vec2, vec2, vec2, vec2);

    force_func get_force;
    func_ptr update_body;

    int 
        WIDTH = 1366, 
        HEIGHT = 768,
        NUM_THREADS = 16,
        obj = 10000,
        max_frames = 5000,
        graph_height = HEIGHT*3/4
    ;
    
    constexpr double  
        G = 6e-7,
        c = 299'792'458,
        spin = .9,
        stV = 0,
        c_squared = c*c
    ;

    double
        time_passed = 0,
        divergence = 14.f,
        takeover = 100,
        e_squared = .086f,
        dt = .001,
        max_temp = 1
    ;

    float 
        initial_a = 1,

        H_0 = 10e-11,
        omega_m = 0.3111,
        omega_l = 0.6889,
        vacuum_density = 3.35, // gev/m^3
        root_omega = sqrt(vacuum_density)
    ;

    float current_a = 1;


    std::mt19937 rng;
    std::uniform_real_distribution<double> angle_rnd(0.0f, 360.0f);
    std::uniform_real_distribution<double> get_rndX(0.0f, WIDTH);
    std::uniform_real_distribution<double> get_rndY(0.0f, HEIGHT); 
    std::uniform_real_distribution<double> normalized_rnd(0.0f, 1.f); 

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