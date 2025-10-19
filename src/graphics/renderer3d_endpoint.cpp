#define GL_SILENCE_DEPRECATION

#include <iomanip>
#include <iostream>
#include <cmath>
#include <GLFW/glfw3.h>
#include <algorithm>

#include <tbb/parallel_sort.h>
using namespace tbb;

#include "Aster/graphics/3d_graphics.h"
#include "Aster/graphics/color_scale.h"
#include "Aster/building-api/logging.h"


namespace Aster{

namespace Renderer{
/**
* @brief bakes a renderer for a specific simulation
* @param s: simulation to render
*/ 
Renderer::Renderer3d* render(Simulation* s){
    return new Renderer::Renderer3d(s);
}

#define NUM_SEGMENTS 15
#define depth_factor 1.2

//===---------------------------------------------------------===//
// 3d rendering impl                                             //
//===---------------------------------------------------------===//
/**
* @brief shows the axis
*/
Renderer3d* Renderer3d::show_axis(){
    show_axis_b = true;
    return this;
}

/**
* @brief returns true if show_axis has been called
* @returns if show_axis has been called
*/
bool Renderer3d::does_show_axis(){
    return show_axis_b;
}


/**
* @brief draws the axis on screen
*/
void Renderer3d::draw_axis(bool is_back){
    // gets the needed points in simulation space
    static vec3 cube_vertecies[8] = {
        {0,0,0}, {0,1,0}, {1,1,0}, {1,0,0},
        {0,0,1}, {0,1,1}, {1,1,1}, {1,0,1} 
    };

    float aspect_ratio = float(current_width) / current_height;

    vec3 mapped_vs[8];
    vec3 excluded = {0,0,100};

    // maps all the points and then finds the fartherst one from the screen
    // in order to exclude the three verticies correlated to that one 
    // to draw them first since they are behind the bodies and the otehr 9 verticies
    #pragma unroll
    for (int i = 0; i < 8; i++) {
        mapped_vs[i] = map_point(cube_vertecies[i]);
        if (mapped_vs[i].z < excluded.z)
            excluded = mapped_vs[i];
    }
    
    auto connect = [excluded, is_back, mapped_vs, aspect_ratio](int x, int y){
        //exstracts the vectors from the array
        vec3 a = mapped_vs[x], b = mapped_vs[y];

        // check if they should be drawn using isback
        if (!((a == excluded || b == excluded) == is_back)) return;

        //sets the color proportional to the distance
        glColor3f((a.z/2 + b.z/2 +2)/2, 0.0, 0.0);

        // connects the points
        glVertex2f(a.x, a.y * aspect_ratio);
        glVertex2f(b.x, b.y * aspect_ratio);
    };
    
    // sets the line thickness
    glLineWidth(2.0f);

    // begines drawing the segments
    glBegin(GL_LINES);
        connect(0, 1);
        connect(0, 3);
        connect(0, 4);
        connect(1, 5);
        connect(1, 2);
        connect(2, 6);
        connect(2, 3);
        connect(3, 7);
        connect(4, 5);
        connect(5, 6);
        connect(4, 7);
        connect(6, 7);
    glEnd();

}

/**
* @brief maps a point from simulation space onto screen space coordinates
* @param v: vector to map
* @returns a mapped a vector
*/
vec3 Renderer3d::map_point(vec3 v, bool fixed) {
    // offsets the vector on the center
    v = v* distance;
    v = v - vec3(.5, .5, .5) * distance;
    v = v*gui_scale;

    v = v*2/gui_scale;

    if (!fixed && !(v.x < 1 && v.x > -1 && v.y < float(current_width)/current_height && v.y > -float(current_width)/current_height))
        return {10,10,10};
        
     v = v*gui_scale * gui_scale/2;

    // applies the z-x rotation
    float x1 =  v.x * cos_x_theta + v.z * sin_x_theta;
    float z1 = -v.x * sin_x_theta + v.z * cos_x_theta;

    // applies the z-y rotation
    float y1 = v.y *  cos_y_theta + z1 * sin_y_theta;
    float z2 = -v.y * sin_y_theta + z1 * cos_y_theta;

    // calculates the parallax strenght
    float scale_z = 1 / (distance + z2/ _s -> get_depth()  + _s -> softening); 
    float scale_x = 1 / (distance + x1/ _s -> get_width()  + _s -> softening);
    float scale_y = 1 / (distance + y1/ _s -> get_height() + _s -> softening);

    v = {x1, y1, z2};
    
    if (fixed)
        return v / distance;
    
    return v;
}


/**
* @brief returns if the vector is container inside the simulation
* @param v: vector to check
* @returns wheter if it is in the simulation
*/
bool Renderer3d::is_unitary_bound(vec3 v, vec3 a){
    return v.x >= -a.x && v.x <= a.x && 
           v.y >= -a.y && v.y <= a.y &&
           v.z >= -a.z && v.z <= a.z;
}

void Renderer3d::setup(){
    render3d = &Renderer3d::draw_detailed3d;

    if (this -> _s -> bodies.positions.size() >= 50)
        render3d = &Renderer3d::draw_minimal3d;

    // defines the center of rotation
    rot_center = {
        _s -> get_width() / 2,
        _s -> get_height()/ 2,
        _s -> get_depth() / 2
    };

    // checks if glfw can initialize correctly
    if (critical_if(!glfwInit(), "failed to load glfw")) 
        exit(-1);
 
    // generates the window
    window = glfwCreateWindow(
        1366, 
        768, 
        "Aster's simulation", 
        nullptr, nullptr
    );
    
    // checks if the window has been created 
    if (critical_if(!window, "glfw failed to create a window")) {
        glfwTerminate();
        exit(0);
    }
}

Renderer3d::Renderer3d(Simulation* _s) : _s(_s){
    
    if (critical_if(!_s, "simulation to render is a nullptr"))
        exit(-1);

    this -> setup();
}

/**
* @brief tells the renderer to shwo the window and render the simulation
*/
void Renderer3d::show(){
    if (critical_if(!_s, "simulation to render is a nullptr"))
        exit(-1);
    
    if (!(window))
        this -> setup();

    // changes the context and defines callback related variables
    glfwMakeContextCurrent(window);
    glfwSetWindowUserPointer(window, this);

    // sets up callbacks
	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback); // window rescaling
    glfwSetKeyCallback(window, handle_keyboard_input);  // keyboard inputs
    glfwSetScrollCallback(window, handle_mouse_scroll); // mouse scroll

    // loads the simulation if it hasn't yet
    if (!_s -> has_loaded_yet())
        _s -> load();

    // starts the main loop 
    while (!glfwWindowShouldClose(window)) {
        // clears the screen
        glClear(GL_COLOR_BUFFER_BIT);
        
        
        // if it is not paused it steps the simulation
        if (!paused) 
            body_update_func();

        // shows the axis if it should
        draw_axis(true);
        
        // calls the rendering function
        (this ->*render3d)();

        draw_axis(false);

        // handles mouse events
        handle_displacement();

        // swaps buffers and pulls events
        glfwSwapBuffers(window); 
        glfwPollEvents();
    }

    // after it has done it closes everything
    glfwDestroyWindow(window);
    glfwTerminate(); 
}

/**
* @brief callback function for window resizing
*/
void Renderer3d::framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    // gets the window pointer
	void* ptr = glfwGetWindowUserPointer(window);
    auto* renderer = static_cast<Renderer3d*>(ptr); // converts it into a renderer

    // resets the size
    renderer -> current_width = width;
    renderer -> current_height = height;
    
    // calibrates the viewport
    glViewport(0, 0, width, height);
}

/**
* @brief handles mouse scrolling
*/
void Renderer3d::handle_mouse_scroll(GLFWwindow* window, double xoffset, double yoffset) {
    // gets the window pointer
    void* ptr = glfwGetWindowUserPointer(window);
    auto* renderer = static_cast<Renderer3d*>(ptr); // converts it into a renderer

    if (critical_if(!renderer, "undefined reference to window")) 
        return; 

    // recalculates and clamps the parallax distance
    renderer -> distance += renderer -> distance * static_cast<float>(yoffset) / 200;
}

/**
*  @brief resets the mouse position
*/
void Renderer3d::reset_mouse(){
    clicked = false;
    mouse_init_x = 0;
    mouse_init_y = 0;
}

/**
* @brief handles the mouse being clicked
*/
void Renderer3d::mouse_clicked(){
    if (clicked) return;
    clicked = true;

    // saves the new cursor position
    glfwGetCursorPos(window, &mouse_init_x, &mouse_init_y);
    init_x_theta = x_theta;
    init_y_theta = y_theta;
}

/**
* @brief handles the mouse drag 
*/
void Renderer3d::handle_displacement(){
    // updates the input dictionary
    if (inputs["up"])
        y_theta += DIS_SCALE;
    if (inputs["down"])    
        y_theta -= DIS_SCALE;
    if (inputs["left"])
        x_theta -= DIS_SCALE;
    if (inputs["right"])
        x_theta += DIS_SCALE;

    // gathers the mouse state
    if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS){
        mouse_clicked();
        double mouse_x, mouse_y;

        // updates the mouse position and angles
        glfwGetCursorPos(window, &mouse_x, &mouse_y);       
        x_theta = init_x_theta - atan2((mouse_init_x - mouse_x) / current_width , 1) * 2.5;
        y_theta = init_y_theta + atan2((mouse_init_y - mouse_y) / current_height, 1) * 2.5; 
    } else
    
    // otherwise it resets the mouse
    reset_mouse();

    // precomputes the cosines and sines of the associated angles 
    cos_x_theta = cos(x_theta);
    cos_y_theta = cos(y_theta);

    sin_x_theta = sin(x_theta);
    sin_y_theta = sin(y_theta); 
}

/**
* @brief rendering method to draw thermal simulations
*/
void Renderer3d::draw_termal3d(){
    // sets up glfw configuration
    glBegin(GL_POINTS);
    glPointSize(10);    

    
    vec3 temp = {0,0,0}; // only makes one instance of a vector
    int index;// keeps track of the color index
    for (int i = 0; i < _s -> bodies.positions.size(); ++i){
        // maps the position to 2d space
        temp = map_point(_s -> bodies.get_position_of(i));

        // calculates the associated color
        index = int(get_coloring_index(_s -> bodies.get_temp_of(i))*255);
        glColor3f(color_scale[index][0], color_scale[index][1], color_scale[index][2]);
    
        // draws the point on screen if it is in the simulation
        if (is_unitary_bound(temp))
            glVertex2f( 
                temp.x, 
                temp.y
            );
    }

    glEnd();
}

/*
* @brief fast method to draw every body in the simulation 
*/
void Renderer3d::draw_minimal3d(){

    // sets up the glfw configuration
    glBegin(GL_POINTS);
    glPointSize(10);    

    vec3 temp = {0,0,0}; // makes only one instance of the position vector

    for (int i = 0; i < _s -> bodies.positions.size(); ++i){
        glColor3f(1.0,1.0,1.0);
        // maps the point in 2d space 
        vec3 v = _s -> bodies.positions[i];
        v.x = (v.x + _s -> get_width()) / (2*_s -> get_width()) ;
        v.y = (v.y + _s -> get_height())/ (2*_s -> get_height());
        v.z = (v.z + _s -> get_depth()) / (2*_s -> get_depth()) ;
        v = map_point(v, false);
    
        // draws the point if it's in range
        if (is_unitary_bound(v))
            glVertex2f( 
                v.x, 
                v.y
            );
    }

    glEnd();
}
/**
* @brief slow method to draw every body in the simulation, useful for few body problems
*/
void Renderer3d::draw_detailed3d() { 
    vec3 mapped_pos = {0, 0,0}; // only one instance
    float aspect_ratio = float(current_width) / float(current_height);
    size_t N = _s -> bodies.positions.size();

    std::vector<std::pair<vec3, int>> sorted_pos(N);

    for (int i = 0; i< N; ++i){
        vec3 v = _s -> bodies.positions[i];
        v.x = (v.x + _s -> get_width()) / (2*_s -> get_width()) ;
        v.y = (v.y + _s -> get_height())/ (2*_s -> get_height());
        v.z = (v.z + _s -> get_depth()) / (2*_s -> get_depth()) ;
        
        sorted_pos[i].first = map_point(v, false);
        sorted_pos[i].second = i;
    }

    parallel_sort(sorted_pos.begin(), sorted_pos.end(), [](const std::pair<vec3, int>& a, const std::pair<vec3, int>& b){return a.first.z > b.first.z;});

    for (int i = 0; i < N; ++i) {
        vec3 v = sorted_pos[i].first;

        // fetches the color from rng colors array
        glColor3f(rng_colors[sorted_pos[i].second % 15][0], 
                  rng_colors[sorted_pos[i].second % 15][1], 
                  rng_colors[sorted_pos[i].second % 15][2]); 

        // calculates the radius based on the log_{10} of the mass by some constant
        REAL radius = std::log10(_s -> bodies.get_mass_of(sorted_pos[i].second)) / 200 * gui_scale / (v.z+2)/2 ;

        // draws a circle around the body using a triangle fan
        glBegin(GL_TRIANGLE_FAN);
        // draws the center of the fan
        //glVertex2f(v.x, v.y * aspect_ratio);

        for (int j = 0; j <= NUM_SEGMENTS; j++) {
            // calculates the points's position with polar coordinates where r = radius
            float angle = 2.0f * M_PI * j / NUM_SEGMENTS;
            float vx = (v.x + cos(angle) * radius);
            float vy = (v.y + sin(angle) * radius * aspect_ratio);

            // draws the vertex
            glVertex2f(vx, vy);
        }
        
        // ends the batch
        glEnd(); 
    }
}

/**
*  @brief internal function to handle keyboard inputs 
*/
void Renderer3d::handle_keyboard_input(GLFWwindow* window, int key, int scancode, int action, int mods) {
    void* ptr = glfwGetWindowUserPointer(window);
    auto* renderer = static_cast<Renderer3d*>(ptr);

    if (!renderer)
        throw std::runtime_error("undefined reference to window");

    if (action == GLFW_PRESS) {
        switch (key) {
            case GLFW_KEY_UP:
                renderer -> inputs["up"] = true;
                break;
            case GLFW_KEY_DOWN:
                renderer -> inputs["down"] = true;
                break;
            case GLFW_KEY_LEFT:
                renderer -> inputs["left"] = true;
                break;
            case GLFW_KEY_RIGHT:
                renderer -> inputs["right"] = true;
                break;
            case GLFW_KEY_SPACE:
                renderer -> inputs["space"] = true;
                renderer -> paused = !renderer -> paused;
                break;
            case GLFW_KEY_KP_ADD:
                renderer -> gui_scale += .1;
                break;
            case GLFW_KEY_KP_SUBTRACT:
                renderer -> gui_scale -= .1;
                break;
            default:
                break;
        }
    } else if (action == GLFW_RELEASE) {
        switch (key) {
            case GLFW_KEY_UP:
                renderer -> inputs["up"] = false;
                break;
            case GLFW_KEY_DOWN:
                renderer -> inputs["down"] = false;
                break;
            case GLFW_KEY_LEFT:
                renderer -> inputs["left"] = false;
                break;
            case GLFW_KEY_RIGHT:
                renderer -> inputs["right"] = false;
                break;
            case GLFW_KEY_SPACE:
                renderer -> inputs["space"] = false;
                break;
            default:
                break;
        }
    }
}

/**
* @brief wrapper function to step the simulation
*/
void Renderer3d::body_update_func(){
    _s -> step();
}

}
}