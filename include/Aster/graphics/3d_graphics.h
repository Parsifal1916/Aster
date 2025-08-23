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

    REAL 
        x_theta = 0, 
        y_theta = 0,
        sin_x_theta = 0,
        cos_x_theta = 0,
        sin_y_theta = 0,
        cos_y_theta = 0,
        distance = 10
    ;   

    Renderer3d(Simulation<vec3>* _s);

    /**
    * @brief wrapper function to step the simulation
    */
    void body_update_func();

    /**
    * @brief shows the axis
    */
    Renderer3d* show_axis();

    /**
    * @brief returns true if show_axis has been called
    * @returns if show_axis has been called
    */
    bool does_show_axis();
 
    /*
    * @brief fast method to draw every body in the simulation 
    */
    void draw_minimal3d();

    /**
    * @brief costly method to draw every body in the simulation, useful for few body problems
    */
    void draw_detailed3d();

    /**
    * @brief rendering method to draw thermal simulations
    */
    void draw_termal3d();
 
    /**
    * @brief tells the renderer to shwo the window and render the simulation
    */
    void show();
    
    std::vector<render_func3d> render_modes3d = {
        &Renderer3d::draw_detailed3d,    
        &Renderer3d::draw_minimal3d, 
        &Renderer3d::draw_minimal3d, 
        &Renderer3d::draw_minimal3d, 
        &Renderer3d::draw_termal3d,
        &Renderer3d::draw_minimal3d
    };
    
    // keyboard to input mapping
    std::unordered_map<std::string,  bool> inputs = {
        {"up", false},
        {"down", false},
        {"left", false},
        {"right", false},
        {"space", false}
    };

    private:
    // should it show the axis on screen?
    bool show_axis_b = false;

    // window current size
    int current_width, current_height;

    // mouse position before clicking
    double mouse_init_x = 0, mouse_init_y = 0;

    /**
    * @brief draws the axis on screen
    */
    void draw_axis();

    bool 
        paused = false ,  // has the simulation been paused?
        clicked = false // has the mouse been clicked?
    ;

    /**
    *  @brief resets the mouse position
    */
    void reset_mouse();


    /**
    * @brief handles the mouse being clicked
    */
    void mouse_clicked();

    /**
    * @brief maps a point from simulation space onto screen space coordinates
    * @param v: vector to map
    * @returns a mapped a vector
    */
    vec3 map_point(vec3 v);

    /**
    * @brief returns if the vector is container inside the simulation
    * @param v: vector to check
    * @returns wheter if it is in the simulation
    */
    bool is_unitary_bound(vec3 v);

    /**
    * @brief handles the mouse drag 
    */
    void handle_displacement();

    /**
    * @brief callback function for window resizing
    */
    static void framebuffer_size_callback(GLFWwindow* window, int width, int height);
    
    /**
    *  @brief internal function to handle keyboard inputs 
    */
    static void handle_keyboard_input(GLFWwindow* window, int key, int scancode, int action, int mods);

    /**
    * @brief handles mouse scrolling
    */
    static void handle_mouse_scroll(GLFWwindow* window, double xoffset, double yoffset);
};




//    void clear_graph3d(Simulation3d* _s);
//    void draw_graph3d(Simulation3d* _s);
}
Renderer::Renderer3d* render(Simulation<vec3>*);

}
