#pragma once
#include <cassert>
#include <chrono>
#include <vector>
#include <cmath>
#include <unordered_map>
#include <functional>
#include <string>

#include "Aster/simulations/3d_sim_obj.h"
    
#define DIS_SCALE .05
#define FOV 1

namespace Aster{
namespace Renderer{
    


class Renderer3d{
    public: 
    bool paused = false;
    using render_func3d = void(Renderer3d::*)();
    Simulation3d* _s = nullptr;
    render_func3d render3d = nullptr;

    GLFWwindow* window;

    vec3 rot_center = vec3(0,0,0); 

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

    Renderer3d(Simulation3d* _s);
    void body_update_func();

 
    void draw_minimal3d();
    void draw_detailed3d();
    void draw_termal3d();
    void run();
    
    std::vector<render_func3d> render_modes3d = {
        &Renderer3d::draw_minimal3d, 
        &Renderer3d::draw_minimal3d, 
        &Renderer3d::draw_minimal3d, 
        &Renderer3d::draw_minimal3d, 
        &Renderer3d::draw_termal3d
    };
    
    std::unordered_map<std::string,  bool> inputs = {
        {"up", false},
        {"down", false},
        {"left", false},
        {"right", false},
        {"space", false}
    };

    private:

    vec3 rotate_point(vec3 v);
    bool is_in3d_bounds(vec3 v);
    void handle_displacement();
    static void handle_keyboard_input(GLFWwindow* window, int key, int scancode, int action, int mods);
    static void handle_mouse_scroll(GLFWwindow* window, double xoffset, double yoffset);
};




//    void clear_graph3d(Simulation3d* _s);
//    void draw_graph3d(Simulation3d* _s);
}
}
