#pragma once
#include <cassert>
#include <SFML/Graphics.hpp>
#include <vector>
#include <unordered_map>
#include <iomanip>

#include "Aster/graphics/clienting.h"
#include "Aster/graphics/inferno_scale.h"
#include "Aster/simulations/sim_obj.h"


namespace Aster{
namespace Renderer{

    extern std::unordered_map<Simulation*, Client> clients;
    
    /*
    * draws a circle of radius $radius and center = (centerX, centerY)
    * note: it draws it on the "image" object
    * @param radius: radius of the circle
    * @param centerX: x component
    * @param centerY: Y component
    */
    void draw_circle(int radius, vec2 center, Simulation* _s, sf::Color color = sf::Color::White);
    void clear_graph(Simulation* _s);
    void draw_graph(Simulation* _s);
    void reset(Simulation* _s);
    int clamp_rgb(double x);

    /*
    * draws every body in the Simulation::bodies
    * with the respective color on the "image"
    * object
    */
    void draw_termal(Simulation* _s);

    /*
    * draws every body in the Simulation::bodies
    * with white on the "image"
    * object
    */
    void draw_minimal(Simulation* _s);
    void draw_detailed(Simulation* _s);
    void serve_new(Simulation* _s);
    inline static bool is_open(Simulation* _s);

    void do_window(Simulation* _s);
        //render_txture.clear();
        //render_txture.draw(sprite);
        //render_txture.draw(Renderer::lagr_text);
        //render_txture.display();
    void save_stepbystep(Simulation* _s);

}
}
