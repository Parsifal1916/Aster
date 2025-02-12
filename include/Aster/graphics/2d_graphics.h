#pragma once
#include <cassert>
#include <vector>
#include <unordered_map>
#include <iomanip>
#include <GLFW/glfw3.h>
#include "Aster/simulations/sim_obj.h"

namespace Aster{

namespace Renderer{

extern float rng_colors[15][3];

class Renderer2d{
    public: 
    using render_func = void(Renderer2d::*)();

    GLFWwindow* window;
    int current_height, current_width;

    Renderer2d(Simulation<vec2>* _s);

    Renderer2d* show_axis();
    bool does_show_axis();

    void body_update_func();
    void draw_termal();
    void draw_minimal();
    void draw_detailed();
    void draw_composite();
    void show();
    static void framebuffer_size_callback(GLFWwindow* window, int width, int height);
    
    std::vector<render_func> render_modes = {
        &Renderer2d::draw_detailed,      
        &Renderer2d::draw_minimal,
        &Renderer2d::draw_minimal,
        &Renderer2d::draw_minimal,
        &Renderer2d::draw_termal,
        &Renderer2d::draw_composite,
    };

    private:

    void draw_axis();

    bool show_axis_b = false;    
    Simulation<vec2>* _s = nullptr;
    render_func render = nullptr;
    
};

}

Renderer::Renderer2d* render(Simulation<vec2>*);

}