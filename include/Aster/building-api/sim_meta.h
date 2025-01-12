#pragma once

#include "Aster/simulations/basic.h"

namespace Aster{

enum simulation_types: int {LIGHT = 0, MEDIUM = 1, HEAVY = 2, BARNES_HUT = 3, BH_termal = 4};


struct sim3d_meta{
   int 
        NUM_THREADS = 16,
        max_frames = 5000
    ;
    
    double  
        G = 6e-11,
        c = 299'792'458,
        c_squared = c*c,
        boltzmann = 5.67e-8,
        pi = 3.141592,
        simulation_scale = 2, // meters / meters
        WIDTH = 1366, HEIGHT = 768,
        depth = 768,
        graph_height = 3/4
    ;

    double
        divergence = 14.f,
        takeover = 100,
        e_squared = 10e-6,
        dt = .001,
        max_temp = 1,
        curvature_factor = 1,
        field_eq_scalar = 8 * pi * G / (c*c*c*c),
        cosmological_constant = 3 * (H_0*H_0*H_0 / (c*c*c))
    ;

    float 
        initial_a = 1,
        current_a = 10,
        H_0 = 10,
        omega_m = 0.3111,
        omega_l = 0.6889,
        vacuum_density = 3.35, // gev/m^3 ? 
        root_omega = sqrt(vacuum_density)
    ;
    bool 
        flrw_on = false,
        save = false,
        show = true,
        leave_traces = false,
        show_graph = false
    ; 
    simulation_types type = HEAVY;

    std::string selected_update = "Leapfrog";
    std::string selected_force = "Pn2"; 
};

struct sim_meta{
   int 
        NUM_THREADS = 16,
        max_frames = 5000
    ;
    
    double  
        G = 6e-11,
        c = 299'792'458,
        c_squared = c*c,
        boltzmann = 5.67e-8,
        pi = 3.141592,
        simulation_scale = 2, // meters / meters
        graph_height = 3/4,
        WIDTH = 1366, 
        HEIGHT = 768
    ;

    double
        divergence = 14.f,
        takeover = 100,
        e_squared = 10e-6,
        dt = .001,
        max_temp = 1,
        curvature_factor = 1,
        field_eq_scalar = 8 * pi * G / (c*c*c*c),
        cosmological_constant = 3 * (H_0*H_0*H_0 / (c*c*c))
    ;

    float 
        initial_a = 1,
        current_a = 10,
        H_0 = 10,
        omega_m = 0.3111,
        omega_l = 0.6889,
        vacuum_density = 3.35, // gev/m^3 ? 
        root_omega = sqrt(vacuum_density)
    ;
    bool 
        flrw_on = false,
        save = false,
        show = true,
        leave_traces = false,
        show_graph = false
    ; 
    simulation_types type = HEAVY;

    std::string selected_update = "Leapfrog";
    std::string selected_force = "Pn2"; 
};

}