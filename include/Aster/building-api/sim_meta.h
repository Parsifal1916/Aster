#pragma once

#include "Aster/simulations/basic.h"

namespace Aster{

enum simulation_types: int {LIGHT = 0, MEDIUM = 1, HEAVY = 2, BARNES_HUT = 3, BH_termal = 4};

struct sim_meta{
   int 
        NUM_THREADS = 16,
        max_frames = 5000
    ;
    
    vec3 size = {
        1366, // x
        768,  // y
        1000  // z
    };

    double  
        avr_heat_capacity = 2e3,
        H_0 = 10,
        G = 6e-11,
        c = 299'792'458,
        c_squared = c*c,
        boltzmann = 1.380649e-23,
        pi = 3.141592,
        simulation_scale = 2, // meters / meters
        graph_height = 3/4
    ;

    double
        divergence = 14.f,
        takeover = 100,
        e_squared = 1,
        dt = .001,
        max_temp = 1,
        curvature_factor = 1,
        field_eq_scalar = 8 * pi * G / (c*c*c*c),
        cosmological_constant = 3 * (H_0*H_0*H_0 / (c*c*c))
    ;

    float 
        initial_a = 1,
        current_a = 10,

        omega_m = 0.3111,
        omega_l = 0.6889,
        vacuum_density = 3.35, // gev/m^3 ? 
        root_omega = sqrt(vacuum_density)
    ;
    bool 
        flrw_on = false,
        save = false
    ; 
    simulation_types type = HEAVY;

    update_type selected_update;
    force_type selected_force; 
};

}