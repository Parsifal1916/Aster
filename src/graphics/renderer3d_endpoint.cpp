#define GL_SILENCE_DEPRECATION

#include <iomanip>
#include <iostream>
#include <cmath>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <algorithm>

#include <tbb/parallel_sort.h>
using namespace tbb;

#include "Aster/graphics/3d_graphics.h"
#include "Aster/graphics/color_scale.h"
#include "Aster/building-api/logging.h"
#include "Aster/graphics/text_rendering.h"

namespace Aster{

namespace Renderer{

#define GRAPH_PADDING .05
#define GRAPH_HEIGHT .3
#define GRAPH_HPAD .4

//===---------------------------------------------------------===//
// SIDEBAR IMPLEMENTATION                                        //
//===---------------------------------------------------------===//
/** 
* @brief draws a rectangle on screen
* @param w width of the rectangle
* @param h height of the rectangle
* @param x x coordinate of the bottom left vertex
* @param y y coordinate of the bottom left vertex
*/
inline void draw_rect(float w, float h, float x, float y){
    glBegin(GL_QUADS);
        glVertex2f( x    , y    ); 
        glVertex2f( x + w, y    ); 
        glVertex2f( x + w, y + h); 
        glVertex2f( x    , y + h); 

    glEnd();    
}

/** 
* @brief A struct that draws graphs from the Graphs::Graph field
* @param s simulation to render
* @param g assigned graph pointer
* @param _r parent renderer
* @param p slot number in the rendering 
*/
GraphDrawer::GraphDrawer(Simulation* s, Graphs::Graph* g, Renderer3d* _r, int p) 
    : _s(s), rend(_r), place(p), _g(g){
    // sets up the default settings for the cube
    corner = {GRAPH_PADDING, GRAPH_HPAD - .4*place};
    data.resize(s -> N);

    // preresizes all the data object
    for (auto& a : data) 
        a.resize(max_size); 
}


/**
 * @brief parses the incoming data from the simulation
 * @param which what body's graph it is rendering
 */
void GraphDrawer::parse_data(int which){
    // returns if the given graph is empty
    if (_g->data[which].size() == 0) return;

    // resets the internal counter if it overflows
    if (internal_counter >= 500) 
        internal_counter = 0;
    
    // copies the data from the simulations
    // graphs to the internal data array
    data[which][index] = _g -> data[which][internal_counter];
 
    // checks for highest and lowest data
    if (data[which][index] > highest) 
        highest = data[which][index];
    
    if (data[which][index] < lowest) 
        lowest = data[which][index]; 
}

/**
 * @brief draws the background of the graph
 */
void GraphDrawer::draw_bg(){
    // saves the height and width of the window
    int w = rend->get_w()/2;
    int h = rend->get_h();

    // sets up the graph's background
    glColor3f(GRAPH_COLOR); 
    draw_rect(1 - 2*GRAPH_PADDING, GRAPH_HEIGHT, corner.x, corner.y);

    // functions to convert form [-1,1] coordinates to pixels
    auto topanel_x = [w](float x){return w * (1 + x);};
    auto topanel_y = [h](float y){return h * (y);};

    //  if no graph is assigned to the renderer it skips the drawing part
    if (_g == nullptr) {
        //  assigns a placeholder text to indicate it's empty as such
        //? "Graphs slot #slot is Empty"
        std::string empty_slot = "Graph Slot " + std::to_string(place+1) + std::string(" is Empty");

        // finally writes the text to screen 
        Text::write2screen(topanel_x(.35), topanel_y(.23 + .2*place), 15, empty_slot.c_str());
        return;
    } 

    // writes the title text of the graph with a custom name
    std::string txt = std::string("Slot ") + std::string(" ") + _g -> name;
    Text::write2screen(w * (1.01+GRAPH_PADDING),  h * (.14 + .2*place) , 15, txt.c_str());
    glColor3f(0.,0.,1.);    
}


/**
 * @brief rendering function for the graph renderer 
 */
void GraphDrawer::draw(){
    // saves the window's height and lenght
    int w = rend->get_w()/2;
    int h = rend->get_h();
    
    // draws the background window
    draw_bg();

    if (!rend -> paused){
        // if no graph is assigned it skips the drawing part
        if (!this -> _g) return;
        
        // if the data is overflowing it halfs it 
        if (index >= max_size){
            for (int i = 0; i < _g->data.size(); i++){
                for (int j = 0; j < max_size / 2; j++){
                    data[i][j] = data[i][j * 2];
                }
            }
            skip_size *= 2;
            index = max_size/2;
        }
        
        // only parses the data every set amount of steps
        if (!(internal_counter % skip_size)){
            for (int i = 0; i < _g->data.size(); ++i)
                parse_data(i);
            index++;
        }

        // increases the internal counter in the 
        // simulation's graphs regardless of the state 
        // of the skip_size
        internal_counter++;
    }
    
    // computes scale and step size in the rendering phase
    REAL scale = .3;
    float step = (1 - 2 * GRAPH_PADDING) / index; 

    auto rescale_data_point = [this](REAL datapoint) -> REAL {
        return (std::abs(this -> highest - this -> lowest) < 1e-9) 
        ? 1
        : (datapoint - this -> lowest) / (this -> highest - this -> lowest);
    };

    // iterates every sub-graph in the assigned _g
    for (int j = 0; j < _g->data.size(); ++j){
        // sets up color and starting postion of the graph
        REAL cx = corner.x;
        REAL offsetY = corner.y;
        auto color = rng_colors[(j + place) % 15];

        // inits opencl api and sets the color
        glBegin(GL_LINE_STRIP);
        glColor3f(color[0], color[1], color[2]);
        
        // goes through every data point
        for (int i = 0; i < index; ++i) {
            REAL y = offsetY + rescale_data_point(data[j][i]) * scale;  
            glVertex2f(cx, y);
            cx += step; 
        }
        glEnd();
    }

    
}

Sidebar::Sidebar(GLFWwindow* _w, Simulation* s, Renderer3d* rend) : window(_w), _s(s) {
    #pragma unroll
    for (int i = 0; i < 3; i++){
        Graphs::Graph* _g = nullptr;
        if (i < _s->graphs.size()+_s->between_graphs.size()) {
            _g = (i < _s->graphs.size()) 
            ? &_s->graphs[i]
            : &_s->between_graphs[i]
            ;
        } 


        g_drawers[i] = {this -> _s, _g, rend, i};
    }
}

void Sidebar::draw(Renderer3d* rend){
    glColor3f(SB_COLOR); 
    draw_rect(1, 2, 0.0, -1.0);

    #pragma unroll
    for (auto& g : g_drawers)
        g.draw();
}


bool can_be_drawn(vec2 point, std::pair<vec2, vec2> size){
    return !(point.x < size.first.x || point.x > size.second.x
             || point.y < size.first.y || point.y > size.second.y);
}

bool can_be_drawn(vec3 point, std::pair<vec2, vec2> size){
    return !(point.x < size.first.x || point.x > size.second.x
             || point.y < size.first.y || point.y > size.second.y);
}

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

void Renderer3d::add_label(std::string name, label_gen _lbl){
    if (err_if(labels.size() > MAX_LABEL_SIZE, "cannot add more than #10 labels to the same renderer")) return;

    labels.emplace_back(name, _lbl);
}

/**
* @brief returns true if show_axis has been called
* @returns if show_axis has been called
*/
bool Renderer3d::does_show_axis(){
    return show_axis_b;
}
void Renderer3d::draw_all_labels(){
    float w = this -> current_width /2;
    float h = this -> current_height;
    float scale = std::min(current_width / 1366.0, current_height /768.0);
    
    auto parse_string = [this](std::pair<std::string, label_gen> label) -> std::string{
            auto val = (label.second)(this -> _s);
            if (val.size() > MAX_NAME_SIZE)
                val = val.substr(0, MAX_NAME_SIZE);
            if (label.first.size() > MAX_NAME_SIZE)
                return label.first.substr(0, MAX_NAME_SIZE) + std::string(": ") + val;
            
            return label.first + std::string(": ") + val;
    };
    
    switch (layout){
        case SIDEBAR:
        for (int i = 0; i < labels.size(); ++i){
            Text::write2screen(w * (1.01+GRAPH_PADDING),  h * (.73 + i*0.03 *scale), scale * 15, parse_string(labels[i]).c_str());
        }
        break;
        case NORMAL:
            return;
        case PANEL:
            for (int i = 0; i < labels.size(); ++i){
                Text::write2screen(w * (.02),  h * (.03 + 0.02*i * scale), scale * 13, parse_string(labels[i]).c_str());
            }
            break;
        default:
        return;
    }
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
        mapped_vs[i] = map_point(cube_vertecies[i]) + cube_offset;
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
    this -> sidebar = Sidebar(window, _s, this);

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

    auto to_scien = [](float in, bool scien, bool is_int) -> std::string {
        std::ostringstream oss;
        if (scien) oss << std::scientific << std::setprecision(2);
        if (is_int) oss << int(in);
        else oss << in;
        return oss.str();
    };

    add_label("FPS", [this, to_scien](Simulation* _) -> std::string{
        return to_scien(this -> fps, false, true);
    });

    add_label("Cycles/s", [this, to_scien](Simulation* _) -> std::string{
        return to_scien(this -> cycles, false, true);
    });

    add_label("N", [](Simulation* _s) -> std::string{
        return std::to_string(_s -> N);
    });

    add_label("Uses_GPU", [](Simulation* _s) -> std::string{
        return (_s -> uses_GPU()) 
        ? "ON"
        : "OFF";
    });

    add_label("Step size", [to_scien](Simulation* _s) -> std::string{
        return to_scien(_s -> get_dt(), false, false);
    });

    //add_label("")
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

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

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

    if (critical_if(glewInit(), "failed to initilize glew"))
        exit(-1);

    Text::load_font();
    Text::nvg_resize(_s -> get_width(), _s -> get_height());

    // starts the main loop 
    while (!glfwWindowShouldClose(window)) {
        auto start = std::chrono::high_resolution_clock::now();
        if (!paused) 
            body_update_func();
        auto c_end = std::chrono::high_resolution_clock::now();

        // shows the axis if it should
        draw_axis(true);
        Text::begin_vg_frame(this -> current_width, this -> current_height);
        
        // calls the rendering function
        (this ->*render3d)();

        draw_axis(false);

        if (layout == SIDEBAR)
            this -> sidebar.draw(this);
        this -> draw_all_labels(); 

        // handles mouse events
        handle_displacement();

        Text::end_vg_frame(); 

        // swaps buffers and pulls events
        glfwSwapBuffers(window); 
        glfwPollEvents();
        glClear(GL_COLOR_BUFFER_BIT);
        auto end = std::chrono::high_resolution_clock::now();

        if (!((int)_s -> get_time_passed() % 30)){
            fps = 1000/(std::chrono::duration<REAL, std::milli>(end - start)).count();
            cycles = 1000/(std::chrono::duration<REAL, std::milli>(c_end - start)).count();
        }
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
        v = map_point(v, false) + cube_offset;
    
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
        
        sorted_pos[i].first = map_point(v, false) + cube_offset;
        sorted_pos[i].second = i;
    }

    parallel_sort(sorted_pos.begin(), sorted_pos.end(), [](const std::pair<vec3, int>& a, const std::pair<vec3, int>& b){return a.first.z > b.first.z;});

    for (int i = 0; i < N; ++i) {
        vec3 v = sorted_pos[i].first;
        if (!can_be_drawn(v, cube_size)) continue;

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

void Renderer3d::change_layout(){
    auto lint = (static_cast<int>(layout) + 1) % 3;
    this -> layout = static_cast<layout_type>(lint);

    switch (layout){
        case NORMAL:
            cube_offset = {0,0,0};
            cube_size = {{-1,-1}, {1,1}};
            break;
        case SIDEBAR:
            cube_offset = {-.5,0,0};
            cube_size = {{-1,-1}, {0,1}};
            break;
        case PANEL:
            cube_offset = {0,0,0};
            cube_size = {{-1,-1}, {1,1}};
            break;
        default: break;
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
            case GLFW_KEY_L:
                renderer -> change_layout();
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