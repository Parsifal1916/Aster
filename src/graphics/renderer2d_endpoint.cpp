#define GL_SILENCE_DEPRECATION

#include <iomanip>
#include <iostream>
#include <cmath>
#include <GLFW/glfw3.h>

#include "Aster/graphics/3d_graphics.h"
#include "Aster/graphics/2d_graphics.h"
#include "Aster/graphics/color_scale.h"
#include "Aster/graphics/animations.h"

#include "Aster/simulations/BHT_sim.h"
#include "Aster/simulations/barnes-hut.h"

#define NUM_SEGMENTS 15

namespace Aster{
    
Renderer::Renderer2d* render(Simulation<vec2>* s){
    return new Renderer::Renderer2d(s);
}

Renderer::Renderer3d* render(Simulation<vec3>* s){
    return new Renderer::Renderer3d(s);
}

namespace Renderer{

float rng_colors[15][3] = {
    {1.0, 0.0, 0.0},
    {0.0, 1.0, 0.0},
    {0.0, 0.0, 1.0},
    {1.0, 1.0, 0.0},
    {1.0, 0.0, 1.0},
    {1.0, 1.0, 1.0},
    {0.0, 1.0, 1.0},
    {1.0,  .5, 0.0},
    {1.0, 0.0, 0.5},
    {0.5, 1.0, 0.0},
    {0.5, 1.0, 0.5},
    {1.0, 0.5, 1.0},
    {0.5, 1.0, 1.0},
    {0.5, 0.5, 1.0},
    {1.0, 0.5, 0.5},
};

//===---------------------------------------------------------===//
// 2d rendering impl                                             //
//===---------------------------------------------------------===//


using render_func = void(*)(Simulation<vec2>*);


Renderer2d* Renderer2d::show_axis(){
	show_axis_b = true;
	return this;
}

bool Renderer2d::does_show_axis(){
	return show_axis_b;
}

Renderer2d::Renderer2d(Simulation<vec2>* _s)  : _s(_s) {
    render = render_modes[_s -> get_type()];

    if (!glfwInit()) 
        throw std::runtime_error("glfw failed to initialize");


    window = glfwCreateWindow(_s -> get_render_height(), _s -> get_render_width(), "Aster's simulation", nullptr, nullptr);

    
    if (!window) {
        glfwTerminate();
        throw std::runtime_error("glfw failed to create window");
    }
}


void Renderer2d::draw_axis(){
    vec2 origin, point1, point2;
    point1.x = _s -> get_render_height();
    point2.y = _s -> get_render_width();

	glBegin(GL_LINES);
        glColor3f(1.0, 0.0, 0.0);
        glVertex2f(origin.x, origin.y);
        glVertex2f(point1.x, point1.y);

        glColor3f(0.0,1.0,0.0);
        glVertex2f(origin.x, origin.y);
        glVertex2f(point2.x, point2.y);
    glEnd();
}


void Renderer2d::framebuffer_size_callback(GLFWwindow* window, int width, int height) {
	void* ptr = glfwGetWindowUserPointer(window);
    auto* renderer = static_cast<Renderer2d*>(ptr);

    renderer -> current_width = width;
    renderer -> current_height = height;
    glViewport(0, 0, width, height);
}


void Renderer2d::show(){
    glfwMakeContextCurrent(window);
    bool paused = false;

	glfwSetWindowUserPointer(window, this);
	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

	if (!_s -> has_loaded_yet()){
		std::thread t([this](){
			this -> _s -> load();
		});
		
		show_loadingbar(window, _s);
		t.join();
	}
	glfwDestroyWindow(window);

 	window = glfwCreateWindow(_s -> get_render_height(), _s -> get_render_width(), "Aster's simulation", nullptr, nullptr);

    glfwMakeContextCurrent(window);
	glfwSetWindowUserPointer(window, this);
	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    while (!glfwWindowShouldClose(window)) {
        if (!paused) body_update_func(); 

        if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
            paused = !paused;

		if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
            break;


        (this ->* render)();

        glfwSwapBuffers(window); 
        glfwPollEvents();
    }

    glfwDestroyWindow(window);
    glfwTerminate();     
}


void Renderer2d::body_update_func(){
    _s -> step();
}
    

void Renderer2d::draw_termal(){
    glClear(GL_COLOR_BUFFER_BIT);
    glBegin(GL_POINTS);

    glPointSize(10);
    double temp = 0;
    
    for (const auto& p : _s -> bodies){
        temp += p.temp;
        int index = int(get_coloring_index(p.temp)*255);
        glColor3f(color_scale[index][0], color_scale[index][1], color_scale[index][2]);
		glVertex2f(
		    2.f *p.position.x/ (_s -> get_width()) - 1, 
		    2.f * p.position.y/ (_s -> get_height()) - 1
		);


    }

    glEnd();
}


void Renderer2d::draw_minimal(){
    glClear(GL_COLOR_BUFFER_BIT);
    glBegin(GL_POINTS);

    glPointSize(10);
    
    for (const auto& p : _s -> bodies){
        glColor3f(1.0, 1.0, 1.0);
		glVertex2f(
		    2.f * p.position.x/ (_s -> get_width()) - 1, 
		    2.f * p.position.y/ (_s -> get_height()) -1
		);
    }

    glEnd();
}


void Renderer2d::draw_detailed(){        
    glClear(GL_COLOR_BUFFER_BIT);

    glClearColor(0.f, 0.f, 0.f, 1.0f);
    
    for (int i = 0; i < _s -> bodies.size(); i++){
        auto& p = _s -> bodies[i];

        glColor3f(rng_colors[i % 14][0], rng_colors[i % 14][1], rng_colors[i % 14][2]); 
        double radius = std::log10(p.mass)/200;

        vec2 mapped_pos = {
            2.f * p.position.x/ _s -> get_width() - 1, 
		    2.f * p.position.y/ _s -> get_height() - 1
        };

        glBegin(GL_TRIANGLE_FAN);

        glVertex2f(
            mapped_pos.x,
            mapped_pos.y
        );

        for (int j = 0; j <= NUM_SEGMENTS; j++) {
            float angle = 2.0f * M_PI * j / NUM_SEGMENTS;
            float vx = mapped_pos.x + cos(angle) * radius;
            float vy = mapped_pos.y + sin(angle) * radius;
            glVertex2f(vx, vy);
        }
        
       glEnd(); 
    }
}

//===---------------------------------------------------------===//
// text rendering                                                //
//===---------------------------------------------------------===//


/*
void clear_graph(Simulation* _s){
    for (int i = 0; i < _s -> get_height(); i++){
    for (int j = _s -> data.graph_height; j < _s -> get_width(); j++){
        clients[_s].screen.setPixel(i, j, bg_color);
    }
    }
}

void draw_graph(Simulation* _s){
    std::ostringstream scien;
    scien << std::scientific << std::setprecision(3) << _s -> lagrangians.back();
    clients[_s].lagr_text.setString("Lagrangian: "+ scien.str() + "J");
    clients[_s].lagr_text.setPosition(10, _s -> data.graph_height + 10);
    clients[_s].window -> draw(clients[_s].lagr_text);

    for (int i = 0; i < _s -> get_height(); ++i)
        clients[_s].screen.setPixel(i, _s -> data.graph_height, sf::Color::White);

    double p1, p2 = 0;
    int max_height = _s -> get_width() - _s -> data.graph_height;

    for (int index = 1; index < _s -> lagrangians.size(); index++){
        p2 = std::clamp(_s -> get_width() - _s -> lagrangians[index] * max_height / _s -> highest_lagrangian, 0.0, _s -> get_width()-1);
        p1 = std::clamp(_s -> get_width() - _s -> lagrangians[index-1] *  max_height / _s -> highest_lagrangian, 0.0, _s -> get_width()-1);

        //assert(p2 < HEIGHT);
        clients[_s].screen.setPixel(index, p2, sf::Color::White);

        int 
            start = (int)std::min(p1, p2),
            stop = (int)std::max(p1, p2)           
        ;

        for (int i = start; i < stop; i++){
            clients[_s].screen.setPixel(index, i, sf::Color::White);
        }

    }
}

void reset(Simulation* _s){
    clients[_s].screen = clients[_s].blank;
}

int clamp_rgb(double x){
    return (int)std::max(std::min(int(x), 255), 0);
}
*/
}
}