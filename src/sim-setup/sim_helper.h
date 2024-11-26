#pragma once
#include <SFML/Graphics.hpp>
#include "../physics/body.h"
#include <string>
#include <map>
#include <functional>
#include "../graphics-api/canvas.h"

typedef std::string str ; 

namespace Simulation{
    // symplectic 4 constants
    constexpr double c1 = 1 / (2 * (2 - std::pow(2, 1.0/3)));
    constexpr double c2 = (1 - pow(2, 1.0/3)) * c1;
    constexpr double d1 = 1 / (2 - std::pow(2, 1.0/3));
    constexpr double d2 = -std::pow(2, 1.0/3) *d1;   

    size_t body_size = sizeof(Body);
    step_func step;

    // updates a single body's acceleration
    void update_pair(Body* b1){
        for (const auto& b2 : bodies){
            if (&b2 == b1) continue;

            b1 -> acceleration += get_force(
                b1 -> mass, b2.mass,
                b1 -> velocity, b2.velocity,
                b1 -> position, b2.position
            ) / b1 -> mass;
    
        }
        lagrangian += .5 * sqr_magn(b1 -> velocity) * b1 -> mass - magnitude(b1 -> acceleration) * b1 -> mass;
    }

    /*
    * updates using euler
    */
    void update_euler(Body* b){
        b -> velocity += b -> acceleration*dt;
        b -> position += b -> velocity *dt;
        b -> acceleration = {0, 0};
    }

    /*
    * updates using leapfrog
    */
    void update_leapfrog(Body* b){
       b -> position += b -> velocity * dt + b -> acceleration * dt *dt * .5;
       b -> velocity += (b -> acceleration + b -> prev_acc)* dt * .5;
       b -> prev_acc = b -> acceleration;
       b -> acceleration =  {0,0};
    }

    /*
    * updates using symplectic 4
    */
    void update_symplectic4(Body* body){
        body -> position += body -> velocity * c1 *dt;
        

        update_pair(body);
        body -> velocity += body -> acceleration * d1 *dt;
        body -> position += body -> velocity * c2 *dt;
        

        update_pair(body);
        body -> velocity += body -> acceleration * d2 *dt;
        body -> position += body -> velocity * c2 *dt;
        
        update_pair(body);
        body -> velocity += body -> acceleration * d1 *dt;
        body -> position += body -> velocity * c1 *dt;
    }

    /*
    * returns the newtonian apporx. for the body's acceleration
    */
    vec2 newtonian(double m1, double m2, vec2 v1, vec2 v2, vec2 p1, vec2 p2){
        vec2 d = p2 - p1 + vec2(1e-10, 1e-10);
        return normalize(d) *G* m1*m1/ (d.x*d.x+d.y*d.y * magnitude(d));
    }

    std::map<str, func_ptr> update_funcs = {
        {"Euler", update_euler} , 
        {"Leapfrog", update_leapfrog},
        {"Symplectic4", update_symplectic4}
    };

    std::map<str, force_func> force_funcs = {
        {"Newton", newtonian} , 
        {"Pn2", pn2}
    };

    str selected_update = "Leapfrog";
    str selected_force = "Pn2";

    /*
    * it is used by the do_parallel method and updates every
    * body position, velocity and acceleration from index
    * index * obj / NUM_THREADS to index + 1 * obj/ NUM_THREADS
    */
    void update_bundle(unsigned short index){
        unsigned int mult, start, stop;
        mult = obj/NUM_THREADS;
        start = index * mult;
        stop = (index + 1) * mult;
              
        for (int i = start; i < stop; ++i){  
           Body* body = &bodies[i];
           update_pair(body);
           update_body(body);

           if (body -> temp > Simulation::max_temp) Simulation::max_temp = body -> temp;
           body -> temp /= max_temp;
        }
    }

    namespace Renderer{
        using render_func = void(*)();
        render_func render;
    }

    /*
    * initializes the simulation
    */
    void init(){
        lagrangians.reserve(WIDTH / 3);
        threads.reserve(NUM_THREADS);
        
        update_body = update_funcs[selected_update];
        get_force = force_funcs[selected_force];
        
        Renderer::init();
        Renderer::render = Renderer::render_modes[sim_type];
    }

    void use_parallel(){
        for (int i = 0; i < NUM_THREADS; ++i)
            threads.push_back(std::thread(update_bundle, i));

        for (auto& t : threads)
            t.join();

        threads.clear();
    }

    void use_single_thread(){
        for (auto& body : bodies){  
            if(body.still)continue;
            update_pair(&body);
            update_body(&body);

            if (body.temp > Simulation::max_temp) Simulation::max_temp = body.temp;
            body.temp /= max_temp;
        }
    }

    /*
    * updates every object by 1 time step
    */
    void general_step(){
        if (sim_type == LIGHT)
            use_single_thread();
        else 
            use_parallel();
        
        if (lagrangian > highest_lagrangian) highest_lagrangian = lagrangian;
        
        lagrangians.push_back(lagrangian);
        
        if (lagrangians.size() > WIDTH)
            lagrangians.erase(lagrangians.begin());
        
        Simulation::max_temp = 1;
        update_scale();
        lagrangian = 0;
    }
}
