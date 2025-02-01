#include <GLFW/glfw3.h>
#include <thread>
#include <iostream>

#include "Aster/graphics/animations.h"
#include "Aster/simulations/sim_obj.h"

namespace Aster{
namespace Renderer{

void draw_quad(GLFWwindow* window, double percent){

    double 
        bar_width = .6,
        bar_height = .05,
        padding = 0.02
    ;

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
    
    glBegin(GL_QUADS);
        float size = 0.1f; 
        glVertex2f(-bar_width + padding , -bar_height + padding);
        glVertex2f(bar_width* percent + padding - .5, - bar_height + padding);
        glVertex2f(bar_width* percent + padding - .5,   bar_height - padding);
        glVertex2f(-bar_width + padding , bar_height - padding);
    glEnd();  
}

void show_loadingbar(GLFWwindow* window, Simulation<vec2>* _s ){
    while (!_s -> has_loaded_yet()){
        glClear(GL_COLOR_BUFFER_BIT);
            
        if (glfwWindowShouldClose(window)){
            glfwDestroyWindow(window);
            glfwTerminate();
            exit(0);
        }

        if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
            glfwSetWindowShouldClose(window, 1);
        }

        draw_quad(
            window,
            _s -> loading_meta.second
        );

        glfwSwapBuffers(window); 
        glfwPollEvents();
    }
}

}
}