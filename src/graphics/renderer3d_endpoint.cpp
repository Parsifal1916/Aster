#define GL_SILENCE_DEPRECATION

#include <iomanip>
#include <iostream>
#include <cmath>
#include <GLFW/glfw3.h>
#include <algorithm>

#include "Aster/graphics/3d_graphics.h"
#include "Aster/graphics/color_scale.h"


namespace Aster{

namespace Renderer{
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
void Renderer3d::draw_axis(){
    // gets the needed points in simulation space
    vec3 origin = {0, 0, 0};
    vec3 point1 = {_s -> get_width(),      0.0,               0.0        };
    vec3 point2 = {       0.0,      _s -> get_height(),       0.0        };
    vec3 point3 = {       0.0,            0.0,          _s -> get_depth()};

    // maps the point onto the screen
    point1 = map_point(point1); // width
    point2 = map_point(point2); // height
    point3 = map_point(point3); // depth
    origin = map_point(origin); // origin
    
    // sets the line thickness
    glLineWidth(2.0f);

    // begines drawing the segments
    glBegin(GL_LINES);
        // connects the width (x) segment
        glColor3f(1.0, 0.0, 0.0);
        glVertex2f(origin.x, origin.y);
        glVertex2f(point1.x, point1.y);

        // conencts the height (y) segment
        glColor3f(0.0, 1.0, 0.0);
        glVertex2f(origin.x, origin.y);
        glVertex2f(point2.x, point2.y);

        // connects the depth (z) segment
        glColor3f(0.0, 0.0, 1.0);
        glVertex2f(origin.x, origin.y);
        glVertex2f(point3.x, point3.y);

    // ends the bacth
    glEnd();

}

/**
* @brief maps a point from simulation space onto screen space coordinates
* @param v: vector to map
* @returns a mapped a vector
*/
vec3 Renderer3d::map_point(vec3 v) {
    // offsets the vector on the center
    v = v - _s->get_center();

    // applies the z-x rotation
    float x1 = v.x * cos_x_theta + v.z * sin_x_theta;
    float z1 = -v.x * sin_x_theta + v.z * cos_x_theta;

    // applies the z-y rotation
    float y1 = v.y * cos_y_theta + z1 * sin_y_theta;
    float z2 = -v.y * sin_y_theta + z1 * cos_y_theta;

    // calculates the parallax strenght
    float scale = 1 / (distance + z2/ _s -> get_width() + .1);

    // calculates the rescaling factor on x y and z
    float scale_x = _s->get_render_width() / _s->get_width();
    float scale_y = _s->get_render_height() / _s->get_height();

    // maps everything onto [-1, 1]
    float ndc_x = (x1 * scale_x * scale) / (_s->get_render_width() / 2.0f);
    float ndc_y = (y1 * scale_y * scale) / (_s->get_render_height() / 2.0f);

    // returns the transformed vector
    return { ndc_x, ndc_y, z2 };
}

/**
* @brief returns if the vector is container inside the simulation
* @param v: vector to check
* @returns wheter if it is in the simulation
*/
bool Renderer3d::is_unitary_bound(vec3 v){
    return v.x >= -1 && v.x <= 1 && 
           v.y >= -1 && v.y <= 1;
}

Renderer3d::Renderer3d(Simulation<vec3>* _s) : _s(_s){
    if (critical_if(!_s, "simulation to render is a nullptr"))
        exit(-1);

    // sets up the user-defined rendering functions
    render3d = render_modes3d[_s -> get_type()];

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
        _s -> get_render_width(), 
        _s -> get_render_height(), 
        "Aster's simulation", 
        nullptr, nullptr
    );
    
    // checks if the window has been created 
    if (critical_if(!window, "glfw failed to create a window")) {
        glfwTerminate();
        exit(0);
    }

}

/**
* @brief tells the renderer to shwo the window and render the simulation
*/
void Renderer3d::show(){
    if (critical_if(!_s, "simulation to render is a nullptr"))
        exit(-1);

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
        if (show_axis_b)
            draw_axis();
        
        // calls the rendering function
        (this ->*render3d)();

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
    renderer -> distance -= renderer -> distance * static_cast<float>(yoffset) / 1000;
    renderer -> distance = (renderer -> distance < .5) ? .5 : renderer -> distance;
    renderer -> distance = (renderer -> distance > 10 * renderer -> _s -> get_width()) ?  10 * renderer -> _s -> get_width() : renderer -> distance;
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
        x_theta = atan2((mouse_init_x - mouse_x)*6, distance);
        y_theta = atan2(-(mouse_init_y - mouse_y)*6, distance);
    }else 
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
        // maps the point in 2d space 
        temp = map_point(_s -> bodies.get_position_of(i));

        // calculates the color based on parallax
        REAL mult = _s -> get_depth()/std::max((temp.z/8), .3);
        glColor3f(
            mult, mult, mult
        );
    
        // draws the point if it's in range
        if (is_unitary_bound(temp))
            glVertex2f( 
                temp.x, 
                temp.y
            );
    }
 
    glEnd();
}
/**
* @brief costly method to draw every body in the simulation, useful for few body problems
*/
void Renderer3d::draw_detailed3d() { 
    vec3 mapped_pos = {0, 0,0}; // only one instance
    float aspect_ratio = float(current_width) / float(current_height);

    for (int i = 0; i < _s -> bodies.positions.size(); ++i) {
        mapped_pos = map_point(_s -> bodies.get_position_of(i)); //saves the position
         
        // fetches the color from rng colors array
        glColor3f(rng_colors[i % 15][0], 
                  rng_colors[i % 15][1], 
                  rng_colors[i % 15][2]); 

        // calculates the radius based on the log_{10} of the mass by some constant
        REAL radius = std::log10(_s -> bodies.get_mass_of(i)) / 3000;

        // draws a circle around the body using a triangle fan
        glBegin(GL_TRIANGLE_FAN);
        // draws the center of the fan
        glVertex2f(mapped_pos.x / aspect_ratio, mapped_pos.y / aspect_ratio);

        for (int j = 0; j <= NUM_SEGMENTS; j++) {
            // calculates the points's position with polar coordinates where r = radius
            float angle = 2.0f * M_PI * j / NUM_SEGMENTS;
            float vx = (mapped_pos.x + cos(angle) * radius) / aspect_ratio;
            float vy = (mapped_pos.y + sin(angle) * radius) / aspect_ratio;

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