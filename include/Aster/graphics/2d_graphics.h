#pragma once
#include <cassert>
#include <vector>
#include <unordered_map>
#include <iomanip>

#include "Aster/graphics/inferno_scale.h"
#include "Aster/simulations/sim_obj.h"

namespace Aster{
namespace Renderer{
    
class Renderer2d{
    public: 
    using render_func = void(Renderer2d::*)();
    Simulation* _s = nullptr;
    render_func render = nullptr;
    
    GLFWwindow* window;

    Renderer2d(Simulation* _s);

    void body_update_func();
    void draw_termal();
    void draw_minimal();
    void draw_detailed();
    void run();
    
    std::vector<render_func> render_modes = {
        &Renderer2d::draw_minimal, 
        &Renderer2d::draw_minimal, 
        &Renderer2d::draw_minimal, 
        &Renderer2d::draw_minimal, 
        &Renderer2d::draw_termal

    };
};

}
}