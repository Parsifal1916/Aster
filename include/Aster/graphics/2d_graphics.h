#pragma once
#include <cassert>
#include <vector>
#include <unordered_map>
#include <iomanip>
#include "Aster/simulations/sim_obj.h"
#include "Aster/graphs/labeling.h"

namespace Aster{

namespace Renderer{

extern float rng_colors[15][3];

class Renderer2d{
    public: 
    using render_func = void(Renderer2d::*)();

    GLFWwindow* window;
    int current_height, current_width;
    Text::LabelQueue<vec2> labels;

    Renderer2d(Simulation<vec2>* _s);

    /**
    * @brief tells the renderer to show the x and y axis
    * @returns a pointer to the renderer
    */
    Renderer2d* show_axis();

    /**
    * @brief returns true if show axis is true
    * @returns wheter show axis is true
    */
    bool does_show_axis();

    /**
    * @brief function that steps the simulation
    */
    void body_update_func();

    /**
    * @brief rendering function made for thermal simualtions
    */
    void draw_termal();

    /**
    * @brief rendering function made to fit most simulations
    */
    void draw_minimal();

    /**
    * @brief slow rendering function with automatic random coloring and circle drawing
    */
    void draw_detailed();

    /**
    * @brief creates a window and renderes the simulation
    */
    void show();

    /**
    * @brief adds a label to the screen
    * @param gen: generator function for the label
    */
    Renderer2d* add_label(Text::label_gen<vec2> gen);

    /**
    * @brief callback for when the window is resized
    */
    static void framebuffer_size_callback(GLFWwindow* window, int width, int height);
    
    std::vector<render_func> render_modes = {
        &Renderer2d::draw_detailed,      
        &Renderer2d::draw_minimal,
        &Renderer2d::draw_minimal,
        &Renderer2d::draw_minimal,
        &Renderer2d::draw_termal,
    };

    private:

    /**
    * @brief draws the x and y axis, automatically called if show_axis is called
    */
    void draw_axis();

    bool show_axis_b = false;    
    Simulation<vec2>* _s = nullptr;
    render_func render = nullptr;
    
};

}

/**
* @brief bakes a renderer for a specific simulation
* @param s: simulation to render
*/ 
Renderer::Renderer2d* render(Simulation<vec2>*);

}