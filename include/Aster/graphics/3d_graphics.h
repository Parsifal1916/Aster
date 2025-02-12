#pragma once
#include <cassert>
#include <chrono>
#include <vector>
#include <cmath>
#include <GLFW/glfw3.h>
#include <unordered_map>
#include <functional>
#include <string>
    
#include "Aster/physics/vectors.h"
#include "Aster/simulations/sim_obj.h"

#define DIS_SCALE .05
#define FOV 1

namespace Aster{

namespace Renderer{

extern float rng_colors[15][3];

class Renderer3d{
    public: 

    using render_func3d = void(Renderer3d::*)();
    Simulation<vec3>* _s = nullptr;
    render_func3d render3d = nullptr;

    GLFWwindow* window;

    vec3 rot_center = vec3(0,0,0); 

    double 
        x_theta = 0, 
        y_theta = 0,
        sin_x_theta = 0,
        cos_x_theta = 0,
        sin_y_theta = 0,
        cos_y_theta = 0,
        distance = 1000
    ;   

    Renderer3d(Simulation<vec3>* _s);
    void body_update_func();

    Renderer3d* show_axis();
    bool does_show_axis();
 
    void draw_minimal3d();
    void draw_detailed3d();
    void draw_termal3d();
    static void framebuffer_size_callback(GLFWwindow* window, int width, int height);
 
    void show();
    
    std::vector<render_func3d> render_modes3d = {
        &Renderer3d::draw_detailed3d,    
        &Renderer3d::draw_minimal3d, 
        &Renderer3d::draw_minimal3d, 
        &Renderer3d::draw_minimal3d, 
        &Renderer3d::draw_termal3d,
        &Renderer3d::draw_minimal3d
    };
    
    std::unordered_map<std::string,  bool> inputs = {
        {"up", false},
        {"down", false},
        {"left", false},
        {"right", false},
        {"space", false}
    };

    private:
    bool show_axis_b = false;
    int current_width, current_height;
    vec2 mouse_init_pos = {0,0};

    void draw_axis();
    bool paused = false , clicked = false;
    void reset_mouse();
    void mouse_clicked();
    vec3 map_point(vec3 v);
    bool is_unitary_bound(vec3 v);
    void handle_displacement();
    static void handle_keyboard_input(GLFWwindow* window, int key, int scancode, int action, int mods);
    static void handle_mouse_scroll(GLFWwindow* window, double xoffset, double yoffset);
};




//    void clear_graph3d(Simulation3d* _s);
//    void draw_graph3d(Simulation3d* _s);
}
Renderer::Renderer3d* render(Simulation<vec3>*);

}
