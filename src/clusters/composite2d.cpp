#include "Aster/simulations/composite.h"
#include "Aster/physics/composite_body.h"

#include "Aster/graphics/2d_graphics.h"
#include <GLFW/glfw3.h>

namespace Aster{

template <typename T> 
void update_bundle(Simulation<T>* _s, unsigned short index);

Composite2d::Composite2d(){
    this -> data = sim_meta();
    this -> data.type = COMPOSITE2D;
    this -> get_force = get_force_func<vec2>(this -> data.selected_force);
    this -> update_body = get_update_func<vec2>(this -> data.selected_update);
    this -> data.graph_height *= this -> data.size.y;

    this -> threads.reserve(this -> get_cores()); 
    this -> obj = this -> bodies.size();
}

inline void Composite2d::step(){
    for (auto& comp : this -> composites)
        comp.update_bodies(this);
    
    for (auto& body : this -> bodies)
        update_body(&body, this);

    this -> obj = this -> bodies.size(); 
    for (int i = 0; i < this -> get_cores(); ++i)
        this -> threads.push_back(std::thread(update_bundle<vec2>, this, i));

    for (auto& t : threads)
        t.join();

    this -> threads.clear();

    for (auto& graph : this -> between_graphs)
        graph.end_batch();


    this -> trigger_all_graphs();
    this -> time_passed++;
}

namespace Renderer{
    void Renderer2d::draw_composite(){   
        glClear(GL_COLOR_BUFFER_BIT);     
        Composite2d* pointer = reinterpret_cast<Composite2d*>(_s);
        
        for (int i = 0; i < pointer -> composites.size(); i++){
            const auto& comp = pointer -> composites[i];
            glBegin(GL_TRIANGLE_FAN);
    
            glColor3f(rng_colors[i % 14][0], rng_colors[i % 14][1], rng_colors[i % 14][2]); 
    
            glVertex2f(
                2.f * pointer -> bodies[comp.get_start()].position.x/ pointer -> get_width() - 1, 
                2.f * pointer -> bodies[comp.get_start()].position.y/ pointer -> get_height() - 1
            );
    
            for (int j = comp.get_start() + 1; j < comp.get_segments() + comp.get_start() + 1; j++) {
                float vx = 2.f * pointer -> bodies[j].position.x / pointer -> get_width() - 1;
                float vy = 2.f * pointer -> bodies[j].position.y / pointer -> get_height()  - 1;
                glVertex2f(vx, vy);
            }

            if (comp.get_segments() % 2){
                float vx = 2.f * pointer -> bodies[comp.get_start() +1].position.x / pointer -> get_width() - 1;
                float vy = 2.f * pointer -> bodies[comp.get_start() +1].position.y / pointer -> get_height()  - 1;
                glVertex2f(vx, vy);
            }
            
           glEnd(); 
        }
    }
}
    
}

