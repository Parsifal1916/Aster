#pragma once

#include "Aster/simulations/basic.h"
#include "Aster/physics/vectors.h"

namespace Aster{

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

    REAL  
        avr_heat_capacity = 2e3,
        H_0 = 10,
        G = 6e-11,
        c = 299'792'458,
        c_squared = c*c,
        boltzmann = 1.380649e-23,
        pi = 3.1415926535,
        simulation_scale = 1, //[ m ]/ [ m ]
        graph_height = 3/4
    ;

    REAL
        adaptive_coeff = 5,
        divergence = 14,
        takeover = 100,
        e_squared = 1,
        dt = 1,
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

    update_type selected_update;
    force_type selected_force; 
};

}