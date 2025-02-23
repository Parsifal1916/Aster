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

Renderer3d* Renderer3d::show_axis(){
    show_axis_b = true;
    return this;
}

bool Renderer3d::does_show_axis(){
    return show_axis_b;
}



void Renderer3d::draw_axis(){
    vec3 origin = {0, 0, 0};
    vec3 point1 = {_s -> get_width(),      0.0,               0.0        };
    vec3 point2 = {       0.0,      _s -> get_height(),       0.0        };
    vec3 point3 = {       0.0,            0.0,          _s -> get_depth()};

    point1 = map_point(point1);
    point2 = map_point(point2);
    point3 = map_point(point3);
    origin = map_point(origin);
    
    glLineWidth(2.0f);

    glBegin(GL_LINES);
        glColor3f(1.0, 0.0, 0.0);
        glVertex2f(origin.x, origin.y);
        glVertex2f(point1.x, point1.y);

        glColor3f(0.0, 1.0, 0.0);
        glVertex2f(origin.x, origin.y);
        glVertex2f(point2.x, point2.y);

        glColor3f(0.0, 0.0, 1.0);
        glVertex2f(origin.x, origin.y);
        glVertex2f(point3.x, point3.y);
    glEnd();

}

vec3 Renderer3d::map_point(vec3 v) {
    v = v - _s->get_center();

    float x1 = v.x * cos_x_theta + v.z * sin_x_theta;
    float z1 = -v.x * sin_x_theta + v.z * cos_x_theta;

    float y1 = v.y * cos_y_theta + z1 * sin_y_theta;
    float z2 = -v.y * sin_y_theta + z1 * cos_y_theta;

    float scale = 1 / (distance + z2/ _s -> get_width() + .1);


    float scale_x = _s->get_render_width() / _s->get_width();
    float scale_y = _s->get_render_height() / _s->get_height();

    float ndc_x = (x1 * scale_x * scale) / (_s->get_render_width() / 2.0f);
    float ndc_y = (y1 * scale_y * scale) / (_s->get_render_height() / 2.0f);

    return { ndc_x, ndc_y, z2 };
}

bool Renderer3d::is_unitary_bound(vec3 v){
    return v.x >= -1 && v.x <= 1 && 
           v.y >= -1 && v.y <= 1;
}

Renderer3d::Renderer3d(Simulation<vec3>* _s) : _s(_s){
    render3d = render_modes3d[_s -> get_type()];
    rot_center = {
        _s -> get_width() / 2,
        _s -> get_height()/ 2,
        _s -> get_depth() / 2
    };

    if (!glfwInit()) 
        return;

    window = glfwCreateWindow(_s -> get_render_width(), _s -> get_render_height(), "Aster's simulation", nullptr, nullptr);
    
    if (!window) {
        glfwTerminate();
        return;
    }

}


void Renderer3d::show(){
    glfwMakeContextCurrent(window);
    glfwSetWindowUserPointer(window, this);

	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetKeyCallback(window, handle_keyboard_input);
    glfwSetScrollCallback(window, handle_mouse_scroll);

    if (!_s -> has_loaded_yet())
        _s -> load();

    while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT);
        if (!paused) body_update_func();

        if (show_axis_b)
            draw_axis();
        
        (this ->*render3d)();

        handle_displacement();

        glfwSwapBuffers(window); 
        glfwPollEvents();
    }

    glfwDestroyWindow(window);
    glfwTerminate(); 
}

void Renderer3d::framebuffer_size_callback(GLFWwindow* window, int width, int height) {
	void* ptr = glfwGetWindowUserPointer(window);
    auto* renderer = static_cast<Renderer3d*>(ptr);

    renderer -> current_width = width;
    renderer -> current_height = height;
    glViewport(0, 0, width, height);
}

void Renderer3d::handle_mouse_scroll(GLFWwindow* window, double xoffset, double yoffset) {
    void* ptr = glfwGetWindowUserPointer(window);
    auto* renderer = static_cast<Renderer3d*>(ptr);
    
    if (!renderer) throw std::runtime_error("undefined reference to window"); 

    renderer -> distance += renderer -> distance * static_cast<float>(yoffset) / 1000;
    renderer -> distance = (renderer -> distance < .5) ? .5 : renderer -> distance;
    renderer -> distance = (renderer -> distance > 10 * renderer -> _s -> get_width()) ?  10 * renderer -> _s -> get_width() : renderer -> distance;
}

void Renderer3d::reset_mouse(){
    clicked = false;
    mouse_init_pos = {0,0};
}

void Renderer3d::mouse_clicked(){
    if (clicked) return;
    clicked = true;
    glfwGetCursorPos(window, &mouse_init_pos.x, &mouse_init_pos.y);
}

void Renderer3d::handle_displacement(){
    if (inputs["up"])
        y_theta += DIS_SCALE;
    if (inputs["down"])    
        y_theta -= DIS_SCALE;
    if (inputs["left"])
        x_theta -= DIS_SCALE;
    if (inputs["right"])
        x_theta += DIS_SCALE;

    if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS){
        mouse_clicked();
        double mouse_x, mouse_y;

        glfwGetCursorPos(window, &mouse_x, &mouse_y);         
        x_theta = atan2((mouse_init_pos.x - mouse_x)*6, distance);
        y_theta = atan2(-(mouse_init_pos.y - mouse_y)*6, distance);
    }else 
        reset_mouse();



    cos_x_theta = cos(x_theta);
    cos_y_theta = cos(y_theta);

    sin_x_theta = sin(x_theta);
    sin_y_theta = sin(y_theta); 
}


void Renderer3d::draw_termal3d(){
    glBegin(GL_POINTS);

    glPointSize(10);    

    vec3 temp = {0,0,0};
    int index;
    for (const auto& p : _s -> bodies){
        // if it's in the canva range it draws it 
        temp = map_point(p.position);

        index = int(get_coloring_index(p.temp)*255);
        glColor3f(color_scale[index][0], color_scale[index][1], color_scale[index][2]);
    
        if (is_unitary_bound(temp))
            glVertex2f( 
                temp.x, 
                temp.y
            );
    }

    glEnd();
}

/*
* draws every body in the Simulation::bodies
* with white on the "image"
* object
*/
void Renderer3d::draw_minimal3d(){
    glBegin(GL_POINTS);

    glPointSize(10);    

    vec3 temp = {0,0,0};
    for (const auto& p : _s -> bodies){
        // if it's in the canva range it draws it 
        temp = map_point(p.position);
        double mult = _s -> get_render_depth()/std::max(temp.z, .001) + .2;
        glColor3f(
            mult, mult, mult
        );
    
    
        if (is_unitary_bound(temp))
            glVertex2f( 
                temp.x, 
                temp.y
            );
    }

    glEnd();
}

void Renderer3d::draw_detailed3d() { 
    vec3 mapped_pos = {0, 0, 0};

    std::vector<std::pair<int, Body<vec3>>> sorted_bodies;
    for (int i = 0; i < _s->bodies.size(); i++) {
        Body<vec3> b = _s->bodies[i];
        b.position = map_point(b.position);
        sorted_bodies.emplace_back(i, b); 
    }

    std::sort(sorted_bodies.begin(), sorted_bodies.end(), 
              [](const std::pair<int, Body<vec3>>& a, const std::pair<int, Body<vec3>>& b) {
                  return a.second.position.z > b.second.position.z;
              });

    for (auto& [original_index, p] : sorted_bodies) {
        mapped_pos = p.position;

        glColor3f(rng_colors[original_index % 15][0], 
                  rng_colors[original_index % 15][1], 
                  rng_colors[original_index % 15][2]); 

        double radius = std::log10(p.mass) / 3000;

        glBegin(GL_TRIANGLE_FAN);
        glVertex2f(mapped_pos.x, mapped_pos.y);

        for (int j = 0; j <= NUM_SEGMENTS; j++) {
            float angle = 2.0f * M_PI * j / NUM_SEGMENTS;
            float vx = mapped_pos.x + cos(angle) * radius;
            float vy = mapped_pos.y + sin(angle) * radius;
            glVertex2f(vx, vy);
        }
        
        glEnd(); 
    }
}


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

void Renderer3d::body_update_func(){
    _s -> step();
}

}
}