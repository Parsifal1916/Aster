#include <GLFW/glfw3.h>
#include <thread>
#include <iostream>
#define GL_SILENCE_DEPRECATION 

#include "Aster/graphics/animations.h"
#include "Aster/simulations/sim_obj.h"
#include "Aster/building-api/logging.h"

namespace Aster{
namespace Renderer{

/**
* @brief draws a loading bar onto a window with a certain percentage
* @param window: window pointer to draw in
* @param percent: loading bar percentage
*/
void draw_quad(GLFWwindow* window, REAL percent){
    if (critical_if(!window, "can't draw a loading bar onto a nullptr"))
        return;
    if (warn_if(percent < 0, "got unexpected loading value"))
        return;

    REAL 
        bar_width = .6,
        bar_height = .05,
        padding = 0.02
    ;

    // draws the outer white frame
    glBegin(GL_LINES);
        glVertex2f(-bar_width, -bar_height);
        glVertex2f( bar_width, -bar_height);

        glVertex2f( bar_width, -bar_height);
        glVertex2f( bar_width,  bar_height);

        glVertex2f( bar_width, bar_height);
        glVertex2f(-bar_width, bar_height);        

        glVertex2f(-bar_width,  bar_height);  
        glVertex2f(-bar_width, -bar_height);
    glEnd();
    
    // draws the inner fill
    glBegin(GL_QUADS);
        float size = 0.1f; // rescales

        // adds a padding to keep it apart from the frame
        glVertex2f(-bar_width + padding , -bar_height + padding);
        glVertex2f(bar_width* percent + padding - .5, - bar_height + padding);
        glVertex2f(bar_width* percent + padding - .5,   bar_height - padding);
        glVertex2f(-bar_width + padding , bar_height - padding);
    glEnd();  
}

/**
* @brief shows a loading bar while loading a simulation's queue
* @param window: window ptr to load in
* @param _s: simulation to load
*/
void show_loadingbar(GLFWwindow* window, Simulation* _s ){
    if (critical_if(!window, "invalid reference to window"))
        exit(-1);
    
    if (critical_if(!_s, "invalid reference to simulation"))
        return;

    while (!_s -> has_loaded_yet()){ // cycles until is has loaded fully
        glClear(GL_COLOR_BUFFER_BIT); // clear the screen

        // checks for interrupts    
        if (glfwWindowShouldClose(window)){
            glfwDestroyWindow(window);
            glfwTerminate();
            exit(0);
        }

        // checks for the escape button being pressed
        if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
            glfwSetWindowShouldClose(window, 1);
        }

        // draws the loading bar
        draw_quad(
            window,
            _s -> loading_meta.second
        );

        // swaps buffers and pulls events
        glfwSwapBuffers(window); 
        glfwPollEvents();
    }
}

}
}