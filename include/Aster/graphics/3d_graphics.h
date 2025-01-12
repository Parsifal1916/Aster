#pragma once
#include <cassert>
#include <chrono>
#include <vector>
#include <cmath>
#include <unordered_map>
#include <string>
#include <SFML/Graphics.hpp>

#include "Aster/simulations/3d_sim_obj.h"
#include "Aster/graphics/clienting.h"
    
#define DIS_SCALE .05
#define FOV 1

namespace Aster{
namespace Renderer{
    bool paused = false;

    std::unordered_map<Simulation3d*, Client> clients3d;
    std::unordered_map<std::string,  bool> inputs = {
        {"up", false},
        {"down", false},
        {"left", false},
        {"right", false},
        {"space", false}
    };

    double 
        mouse_x = 0, 
        mouse_y = 0, 
        x_theta = 0, 
        y_theta = 0,
        sin_x_theta = 0,
        cos_x_theta = 0,
        sin_y_theta = 0,
        cos_y_theta = 0,
        distance
    ;

    vec3 rot_center = vec3(0,0,0);


    vec3 rotate_point(Simulation3d* _s, vec3 v);
    bool is_in3d_bounds(Simulation3d* _s, vec3 v);
    void clear_graph3d(Simulation3d* _s);
    void draw_graph3d(Simulation3d* _s);

    /*
    * draws every body in the Simulation::bodies
    * with the respective color on the "image"
    * object
    */
    void draw_termal3d(Simulation3d* _s);
    void draw_sphere(int radius, vec3 center, Simulation3d* _s, sf::Color color = sf::Color::White);
    void clear_graph(Simulation3d* _s);

    void handle_displacement();

    void reset(Simulation3d* _s);

    /*
    * draws every body in the Simulation::bodies
    * with white on the "image"
    * object
    */
    void draw_minimal3d(Simulation3d* _s);
    void draw_detailed3d(Simulation3d* _s);

    void serve_new3d(Simulation3d* _s);
    inline static bool is_open(Simulation3d* _s);
    void handle_keyboard_input(sf::Event e);

    void do_window3d(Simulation3d* _s);
        //render_txture.clear();
        //render_txture.draw(sprite);
        //render_txture.draw(Renderer::lagr_text);
        //render_txture.display();
    void save_stepbystep(Simulation3d* _s);

}
}
