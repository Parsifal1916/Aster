#define GL_SILENCE_DEPRECATION

#include <iomanip>
#include <iostream>
#include <cmath>
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "Aster/graphics/2d_graphics.h"
#include "Aster/graphics/color_scale.h"
#include "Aster/graphics/animations.h"
#include "Aster/building-api/logging.h"
#include "Aster/simulations/barnes-hut.h"
#include "Aster/graphics/text_rendering.h"

#define NUM_SEGMENTS 15

namespace Aster{



/**
* @brief bakes a renderer for a specific simulation
* @param s: simulation to render
*/ 
Renderer::Renderer2d* render(Simulation* s){
    return new Renderer::Renderer2d(s);
}



namespace Renderer{

/**
* @brief updates the scale by ns
* @param ns: value to update the scale with
*/
Renderer2d* Renderer2d::update_scale(REAL ns){
    scale += ns;
    return this;
}

/**
* @brief return the scale
*/
REAL Renderer2d::get_scale() const{
    return scale;
}

/**
* @brief handles mouse scrolling
*/
void Renderer2d::handle_mouse_scroll(GLFWwindow* window, double xoffset, double yoffset) {
    // gets the window pointer
    void* ptr = glfwGetWindowUserPointer(window);
    auto* renderer = static_cast<Renderer2d*>(ptr); // converts it into a renderer

    if (critical_if(!renderer, "undefined reference to window")) 
        return; 

    renderer -> update_scale(yoffset/50);
}
/**
* @brief adds a label to the screen
* @param gen: generator function for the label
*/
Renderer2d* Renderer2d::add_label(Text::label_gen gen){
    this -> labels.add_label(gen);
    return this;
}


// used to generate random colors for draw detailed
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


using render_func = void(*)(Simulation*);

/**
* @brief tells the renderer to show the x and y axis
* @returns a pointer to the renderer
*/
Renderer2d* Renderer2d::show_axis(){
	show_axis_b = true;
	return this;
}

/**
* @brief returns true if show axis is true
* @returns wheter show axis is true
*/
bool Renderer2d::does_show_axis(){
	return show_axis_b;
}

Renderer2d::Renderer2d(Simulation* _s)  : _s(_s) {
    if (critical_if(!_s, "simulation to render is a nullptr"))
        exit(-1);
    render = render_modes[0];

    if (critical_if(!glfwInit(), "failed to load glfw")) 
        exit(-1);

    // makes a window based on the simulation's size
    window = glfwCreateWindow(
        _s -> get_render_height(), 
        _s -> get_render_width(), 
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
* @brief draws the x and y axis, automatically called if show_axis is called
*/
void Renderer2d::draw_axis(){
    vec2 origin, point1, point2;

    // creates points on screen representing the max y and x
    point1.x = _s -> get_render_height();
    point2.y = _s -> get_render_width();

    // draws the segments
	glBegin(GL_LINES);
        // draws the x axis with red
        glColor3f(1.0, 0.0, 0.0);
        glVertex2f(origin.x, origin.y);
        glVertex2f(point1.x, point1.y);

        // draws the y axis with green
        glColor3f(0.0,1.0,0.0);
        glVertex2f(origin.x, origin.y);
        glVertex2f(point2.x, point2.y);
    glEnd();
}

/**
* @brief callback for when the window is resized
*/
void Renderer2d::framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    // gets the associated pointer with the window
	void* ptr = glfwGetWindowUserPointer(window);
    // makes it into a renderer pointer
    auto* renderer = static_cast<Renderer2d*>(ptr);

    // resets the renderer width and height
    renderer -> current_width = width;
    renderer -> current_height = height;

    Text::nvg_resize(width, height);

    // sets the right viewport
    glViewport(0, 0, width, height);
}

/**
* @brief creates a window and renderes the simulation
*/
void Renderer2d::show(){
    if (critical_if(!_s, "the given simulation pointer is a nullptr"))
        exit(-1);

    glfwMakeContextCurrent(window);
    bool paused = false;

    // sets up the associated pointer for the callbacks
	glfwSetWindowUserPointer(window, this);
    // stets up the callbacks
	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetScrollCallback(window, handle_mouse_scroll); // mouse scroll

    
    // if the simulation has not been loaded it load it
	if (!_s -> has_loaded_yet()){
		std::thread t([this](){
			this -> _s -> load();
		});
		
        // spawns a loading bar waiting for the thread
		show_loadingbar(window, _s);
		t.join();
	}

    // destroys the loading bar window
	glfwDestroyWindow(window);

    // creates the main rendering window
 	window = glfwCreateWindow(
        _s -> get_render_height(), 
        _s -> get_render_width(), 
        "Aster's simulation", 
        nullptr, nullptr
    );


    // resets the context
    glfwMakeContextCurrent(window);
    // sets up the associated pointer for the callbacks
	glfwSetWindowUserPointer(window, this);
    // stets up the callbacks
	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetScrollCallback(window, handle_mouse_scroll); // mouse scroll

    // initializes glew
    if (critical_if(glewInit(), "failed to initilize glew"))
        exit(-1);

    // initilizes window creation
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    // loads vg and font
    Text::load_font();
    Text::nvg_resize(_s -> get_width(), _s -> get_height());

    // starts the main loop
    while (!glfwWindowShouldClose(window)) {
        // calls the update function of the bodies
        if (!paused) body_update_func(); 

        // checks the space bar for pausing action
        if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
            paused = !paused;

        // if escape is pressed it leaves the simulation
		if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
            break;
        
        // calls the user-defined rendering function
        (this ->* render)();
        
        // loads labels and ends vg frame
        labels.load(_s);
        Text::end_vg_frame(); 

        // swaps buffers and pulls the events 
        glfwSwapBuffers(window); 
        glfwPollEvents();
    }

    // after it has left it destroys everything
    glfwDestroyWindow(window);
    glfwTerminate();     
}

/**
* @brief function that steps the simulation
*/
void Renderer2d::body_update_func(){
    _s -> step();
}
    
/**
* @brief rendering function made for thermal simualtions
*/
void Renderer2d::draw_termal(){
    // clears the screen and starts drawing points
    glClear(GL_COLOR_BUFFER_BIT);
    glBegin(GL_POINTS);

    // sets the point size to 10
    glPointSize(10);
    
    for (int i = 0; i < _s -> bodies.positions.size(); ++i){
        // trasforms the coordinates to -> [-1, 1]
        REAL x = 2.f * (_s -> bodies.get_position_of(i).x * get_scale())/ (_s -> get_width());
        REAL y = 2.f * (_s -> bodies.get_position_of(i).y * get_scale())/ (_s -> get_height());
           
        //in-bounds check
        if (std::abs(x) > 1 || std::abs(y) > 1) continue;

        // gets the right color for the specific temperature
        int index = int(
            get_coloring_index(_s -> bodies.get_temp_of(i))*255
        );
        // changes the color based on the index
        glColor3f(color_scale[index][0], color_scale[index][1], color_scale[index][2]);
		
        // draws the point with the right color
        glVertex2f(x, y);

    }

    // ends the batch
    glEnd();
}


/**
* @brief rendering function made to fit most simulations
*/
void Renderer2d::draw_minimal(){
    // clears the screen and starts drawing points
    glClear(GL_COLOR_BUFFER_BIT);
    glColor3f(1.0, 1.0, 1.0);
    glBegin(GL_POINTS);

    // sets the point size to 10
    glPointSize(10);
    
    for (int i = 0; i < _s -> bodies.positions.size(); ++i){
        // trasforms the coordinates to -> [-1, 1]
        REAL x = 2.f * (_s -> bodies.get_position_of(i).x * get_scale())/ (_s -> get_width());
        REAL y = 2.f * (_s -> bodies.get_position_of(i).y * get_scale())/ (_s -> get_height());
           
        //in-bounds check
        if (std::abs(x) > 1 || std::abs(y) > 1) continue;

        glVertex2f(x, y);
    }

    // ends the batch
    glEnd();
}


/**
* @brief rendering function made for barnes-hut simulations
*/
inline void Renderer2d::draw_barnes(){
    // clears the screen and starts drawing points
    glClear(GL_COLOR_BUFFER_BIT);
    glColor3f(1.0, 1.0, 1.0);
    glBegin(GL_POINTS);

    auto* sim = reinterpret_cast<Barnes::Barnes_Hut*>(_s);

    // sets the point size to 10
    glPointSize(10);

    for (int i = 0; i < _s -> bodies.positions.size(); ++i){
        REAL x = 2.f * (_s -> bodies.get_position_of(i).x * get_scale())/ (_s -> get_width());
        REAL y = 2.f * (_s -> bodies.get_position_of(i).y * get_scale())/ (_s -> get_height());
           
        //in-bounds check
        if (std::abs(x) > 1 || std::abs(y) > 1) continue;

        glVertex2f(x, y);
    }

    glEnd();


}

/**
* @brief slow rendering function with automatic random coloring and circle drawing
*/
void Renderer2d::draw_detailed(){   
    // clears the screen and starts drawing points
    glClear(GL_COLOR_BUFFER_BIT);
    glClearColor(0.f, 0.f, 0.f, 1.0f);
    
    // calculates the current ascpect ratio
    float aspect_ratio = static_cast<float>(current_width) / current_height;
    
    for (int i = 0; i < _s -> bodies.positions.size(); i++){
        // saves a reference to the current body

        // chooses a random color from the array
        glColor3f(rng_colors[i % 14][0], rng_colors[i % 14][1], rng_colors[i % 14][2]); 
        REAL radius = std::max(
            std::min(
                double( std::log10(_s -> bodies.get_mass_of(i) * _s -> get_render_width() / _s -> get_width() ) / (200 * get_scale())),
                .2
            ),
            0.005
        ); // scales the radius based on the log_{10} of the mass

        // trasforms the coordinates to -> [-1, 1]
        REAL x = 2.f * (_s -> bodies.get_position_of(i).x * get_scale())/ (_s -> get_width());
        REAL y = 2.f * (_s -> bodies.get_position_of(i).y * get_scale())/ (_s -> get_height());
           
        //in-bounds check
        if (std::abs(x) > 1.3 || std::abs(y) > 1.3) continue; 

        /*
        * the 1.3 in the if above is there to prevent objects form suddenly disappearing 
        * when almost on the edge of the screen
        */

        // calculates the center 
        vec2 mapped_pos = {x, y};

        // begins drawing the circle
        glBegin(GL_TRIANGLE_FAN);
        // draws the center
        glVertex2f(mapped_pos.x, mapped_pos.y);

        for (int j = 0; j <= NUM_SEGMENTS; j++) {
            // calculates the next segment with polar coordinates
            float angle = 2.0f * M_PI * j / NUM_SEGMENTS;
            // evaluates the x componenet
            float vx = mapped_pos.x + get_scale() * radius * cos(angle) / (aspect_ratio > 1 ? aspect_ratio : 1);
            // evaliates the y compoment
            float vy = mapped_pos.y + get_scale() * radius * sin(angle) * (aspect_ratio < 1 ? aspect_ratio : 1);
            glVertex2f(vx, vy);
        }
        
        // ends the batch
        glEnd(); 
    }
}


//===---------------------------------------------------------===//
// text rendering                                                //
//===---------------------------------------------------------===//


/*
! hot garbage 
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

    REAL p1, p2 = 0;
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

int clamp_rgb(REAL x){
    return (int)std::max(std::min(int(x), 255), 0);
}
*/
}
}