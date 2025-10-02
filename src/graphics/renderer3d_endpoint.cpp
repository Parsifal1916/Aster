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
    vec3 origin = {0.0, 0.0, 0.0};
    vec3 opp =    {1.0, 1.0, 1.0};
    vec3 point1 = {1.0, 0.0, 0.0};
    vec3 point2 = {0.0, 1.0, 0.0};
    vec3 point3 = {0.0, 0.0, 1.0};
    vec3 mid1   = {1.0, 1.0, 0.0};
    vec3 mid2   = {0.0, 1.0, 1.0};
    vec3 mid3   = {1.0, 0.0, 1.0};

    // maps the point onto the screen
    point1 = map_point(point1); // width
    point2 = map_point(point2); // height
    point3 = map_point(point3); // depth
    origin = map_point(origin); // origin
    mid1 = map_point(mid1);
    mid2 = map_point(mid2);
    mid3 = map_point(mid3);
    opp = map_point(opp);

    auto connect = [](vec3 a, vec3 b){
        glColor3f(1.0, 0.0, 0.0);
        glVertex2f(a.x, a.y);
        glVertex2f(b.x, b.y);
    };
    
    // sets the line thickness
    glLineWidth(2.0f);

    // begines drawing the segments
    glBegin(GL_LINES);
        connect(origin, point1);
        connect(origin, point2);
        connect(origin, point3);

        connect(mid2, point2);
        connect(mid2, point3);
        
        connect(mid1, point1);
        connect(mid1, point2);

        connect(mid3, point1);
        connect(mid3, point3);

        connect(mid1, opp);
        connect(mid2, opp);
        connect(mid3, opp);
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

    if (!fixed && !is_unitary_bound(v*2))
        return {10,10,10};
        
    

    // applies the z-x rotation
    float x1 =  v.x * cos_x_theta + v.z * sin_x_theta;
    float z1 = -v.x * sin_x_theta + v.z * cos_x_theta;

    // applies the z-y rotation
    float y1 = v.y *  cos_y_theta + z1 * sin_y_theta;
    float z2 = -v.y * sin_y_theta + z1 * cos_y_theta;

    // calculates the parallax strenght
    float scale_z = 1 / (distance + z2/ _s -> get_depth()  + 1e-10); 
    float scale_x = 1 / (distance + x1/ _s -> get_width()  + 1e-10);
    float scale_y = 1 / (distance + y1/ _s -> get_height() + 1e-10);

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
        x_theta = atan2((mouse_init_x - mouse_x), distance);
        y_theta = atan2(-(mouse_init_y - mouse_y), distance);
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