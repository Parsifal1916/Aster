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
#define MAX_LABEL_SIZE 10
#define MAX_NAME_SIZE 16
#define SB_COLOR 48.0/255.0, 48.0/255.0, 48.0/255.0
#define GRAPH_COLOR 28.0/255.0, 28.0/255.0, 28.0/255.0

namespace Aster{

namespace Renderer{
enum layout_type:  int {NORMAL = 0, SIDEBAR = 1, PANEL = 2};
extern float rng_colors[15][3];
struct Renderer3d;
using label_gen = std::function<std::string(Simulation*)>;

struct GraphDrawer{
    GraphDrawer(Simulation* s, Graphs::Graph* g, Renderer3d* _r, int place);
    GraphDrawer(){};
    void draw();
    void draw_bg();
    void parse_data(int which);

    int index = 0;
    int skip_size = 1;
    vec2 corner;
    REAL highest = std::numeric_limits<double>::lowest();
    REAL lowest = std::numeric_limits<double>::max();
    int place =0;
    int internal_counter = 0;
    int max_size = 1024;
    std::vector<std::vector<REAL>> data;
    Simulation* _s;
    Renderer3d* rend;
    Graphs::Graph* _g;
};
struct Sidebar {
    Sidebar(GLFWwindow* _w, Simulation* s, Renderer3d* rend);
    Sidebar(){}
    void draw(Renderer3d* rend);
    GLFWwindow* window;
    Simulation* _s;

    GraphDrawer g_drawers[3];
};

class Renderer3d{
    public: 
    using render_func3d = void(Renderer3d::*)();
    Simulation* _s = nullptr;
    render_func3d render3d = &Renderer3d::draw_detailed3d;

    GLFWwindow* window;

    vec3 rot_center = vec3(0,0,0); 

    REAL 
        x_theta = 0, 
        y_theta = 0,
        init_x_theta = 0, 
        init_y_theta = 0,
        sin_x_theta = 0,
        cos_x_theta = 0,
        sin_y_theta = 0,
        cos_y_theta = 0,
        distance = 1
    ;   

    Renderer3d(Simulation* _s);

    /**
    * @brief wrapper function to step the simulation
    */
    void body_update_func();

    void change_layout();

    /**
    * @brief shows the axis
    */
    Renderer3d* show_axis();

    void add_label(std::string name, label_gen _lbl);
    void draw_all_labels();

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
    
    std::vector<render_func3d> render_modes3d;
    
    // keyboard to input mapping
    std::unordered_map<std::string,  bool> inputs = {
        {"up", false},
        {"down", false},
        {"left", false},
        {"right", false},
        {"space", false}
    };
    REAL gui_scale = .5;

    int get_h() const {return this -> current_height;}
    int get_w() const {return this -> current_width;}

    bool 
        paused = false ,  // has the simulation been paused?
        clicked = false // has the mouse been clicked?
    ;

    std::vector<std::pair<std::string, label_gen>> labels;
    private:
    float fps = 0.0;
    float cycles  = 0.0;
    Sidebar sidebar;
    vec3 cube_offset = {0,0, 0};
    std::pair<vec2, vec2> cube_size = {{-1,-1}, {1,1}};
    layout_type layout = NORMAL;


    void setup();
    // should it show the axis on screen?
    bool show_axis_b = true;

    // window current size
    int current_width, current_height;

    // mouse position before clicking
    double mouse_init_x = 0, mouse_init_y = 0;

    /**
    * @brief draws the axis on screen
    */
    void draw_axis(bool is_back);

 

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
    vec3 map_point(vec3 v, bool fixed = true);

    /**
    * @brief returns if the vector is container inside the simulation
    * @param v: vector to check
    * @returns wheter if it is in the simulation
    */
    bool is_unitary_bound(vec3 v, vec3 a = {1.0,1.0,1.0});

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
Renderer::Renderer3d* render(Simulation*);

}
