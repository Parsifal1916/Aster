#define GL_SILENCE_DEPRECATION

#include <iomanip>
#include <iostream>
#include <cmath>
#include <GLFW/glfw3.h>

#include "Aster/graphics/3d_graphics.h"
#include "Aster/graphics/inferno_scale.h"

#include "Aster/simulations/barnes-hut3d.h"

namespace Aster{
    
void render(Simulation3d* s){
    auto r = Renderer::Renderer3d(s);
    r.run();
}

namespace Renderer{

    #define depth_factor 1.2

//===---------------------------------------------------------===//
// 3d rendering impl                                             //
//===---------------------------------------------------------===//

vec3 Renderer3d::rotate_point(vec3 v){
    v =  v - _s -> get_center(); 

    v.x /= _s -> data.WIDTH  / 2 - 1;
    v.y /= _s -> data.HEIGHT / 2 - 1;
    v.z /= _s -> data.depth  / 2 - 1;

    float x1, z1, y1, z2;

    x1 = v.x * cos_x_theta + v.z * sin_x_theta;
    z1 = - v.x * sin_x_theta + v.z *cos_x_theta;

    y1 = v.y * cos_y_theta + z1 *sin_y_theta;
    z2 = -v.y *sin_y_theta + z1 * cos_y_theta;


    float scale = distance/ std::pow(distance + z2, depth_factor);

    v.x = (x1*scale) * (double)window_width  / _s -> data.WIDTH ; 
    v.y = (y1*scale) * (double)window_height / _s -> data.HEIGHT; 
    v.z = z2;

    return v;
}

bool Renderer3d::is_unitary_bound(vec3 v){
    return v.x >= -1 && v.x <= 1 && 
           v.y >= -1 && v.y <= 1;
}

Renderer3d::Renderer3d(Simulation3d* _s) : _s(_s){
    render3d = render_modes3d[_s -> data.type];
    rot_center = {
        _s -> data.WIDTH / 4,
        _s -> data.HEIGHT/ 4,
        _s -> data.depth / 4
    };

    if (!glfwInit()) 
        return;

    window = glfwCreateWindow(_s -> data.WIDTH, _s -> data.HEIGHT, "Aster's simulation", nullptr, nullptr);
    
    if (!window) {
        glfwTerminate();
        return;
    }

}


void Renderer3d::run(){
    //glViewport(-_s -> data.WIDTH/2, -_s -> data.HEIGHT/2, -_s -> data.WIDTH, -_s -> data.HEIGHT);
    
    glfwMakeContextCurrent(window);
    glfwSetWindowUserPointer(window, this);

    glfwSetKeyCallback(window, handle_keyboard_input);
    glfwSetScrollCallback(window, handle_mouse_scroll);

    while (!glfwWindowShouldClose(window)) {
        glfwGetFramebufferSize(window, &window_width, &window_height);
        if (!paused) body_update_func();

        (this ->*render3d)();

        handle_displacement();

        glfwSwapBuffers(window); 
        glfwPollEvents();
    }

    glfwDestroyWindow(window);
    glfwTerminate(); 
}

void Renderer3d::handle_mouse_scroll(GLFWwindow* window, double xoffset, double yoffset) {
    void* ptr = glfwGetWindowUserPointer(window);
    auto* renderer = static_cast<Renderer3d*>(ptr);
    
    if (!renderer) throw std::runtime_error("undefined reference to window"); 

    renderer -> distance += static_cast<float>(yoffset) * 100.0f * renderer -> distance / (renderer -> distance +1000);
    renderer -> distance = (renderer -> distance < .5) ? .5 : renderer -> distance;
    renderer -> distance = (renderer -> distance > 10 * renderer -> _s -> data.WIDTH) ?  10 * renderer -> _s -> data.WIDTH : renderer -> distance;
}

void Renderer3d::draw_termal3d(){
    glClearColor(0.f, 0.f, 0.f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    glBegin(GL_POINTS);

    glPointSize(10);    

    vec3 temp = {0,0,0};
    for (const auto& p : _s -> bodies){
        glColor3f(1.0, 1.0, 1.0);
        temp = rotate_point(p.position);
        if (is_unitary_bound(temp))        // if it's in the canva range it draws it 
            glVertex2f(
                temp.x, 
                temp.y
            );
    }

    glEnd();
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

/*
* draws every body in the Simulation::bodies
* with white on the "image"
* objectc
*/
void Renderer3d::draw_minimal3d(){
    glClearColor(0.f, 0.f, 0.f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    glBegin(GL_POINTS);

    glPointSize(10);    

    vec3 temp = {0,0,0};
    for (const auto& p : _s -> bodies){
        glColor3f(1.0, 1.0, 1.0);
        // if it's in the canva range it draws it 
        temp = rotate_point(p.position);
        if (is_unitary_bound(temp))
            glVertex2f(
                temp.x, 
                temp.y
            );
    }

    glEnd();
}

void Renderer3d::draw_detailed3d(){
    glClearColor(0.f, 0.f, 0.f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    glBegin(GL_POINTS);

    vec3 temp = {0,0,0};
    for (const auto& p : _s -> bodies){
        glColor3f(1.0, 1.0, 1.0);
        glPointSize(std::log(p.mass)/10);
        // if it's in the canva range it draws it 
        temp = rotate_point(p.position);
        if (is_unitary_bound(temp))
            glVertex2f(
                temp.x, 
                temp.y
            );
    }

    glEnd();
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
    _s -> time_passed++;
}

}
}