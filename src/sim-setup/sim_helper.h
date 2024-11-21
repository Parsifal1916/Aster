#pragma once
#include <SFML/Graphics.hpp>
#include "../physics/body.h"

namespace Simulation{
    size_t body_size = sizeof(Body);

    void update_pair(struct Body& b1){
        for (const auto& b2 : bodies){
            if (&b2 == &b1) continue;

            b1.acceleration += pn2(
                b1.mass, b2.mass,
                b1.velocity, b2.velocity,
                b1.position, b2.position
            );
        }
    }

    void update_bundle(unsigned short index){
        unsigned int mult, start, stop;
        mult = obj/NUM_THREADS;
        start = index * mult;
        stop = (index + 1) * mult;
              
        for (int i = start; i < stop; ++i){  
           Body& body = bodies[i];
           update_pair(body);

           body.velocity += body.acceleration*dt;
           body.position += body.velocity *dt;
           body.acceleration = {0, 0};

           if (body.temp > Simulation::max_temp) Simulation::max_temp = body.temp;
           body.temp /= max_temp;
           body.color = Body::get_color(body.temp);
        }
    }

    void init(){
        threads.reserve(NUM_THREADS);
    }

    void step(){
        for (int i = 0; i < NUM_THREADS; ++i)
            threads.push_back(std::thread(update_bundle, i));

        for (auto& t : threads)
            t.join();

        threads.clear();
            
        Simulation::max_temp = 1;
        update_scale();
    }
}
