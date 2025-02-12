// *  █████╗ ███████╗████████╗███████╗██████╗ 
// * ██╔══██╗██╔════╝╚══██╔══╝██╔════╝██╔══██╗
// * ███████║███████╗   ██║   █████╗  ██████╔╝
// * ██╔══██║╚════██║   ██║   ██╔══╝  ██╔══██╗
// * ██║  ██║███████║   ██║   ███████╗██║  ██║
// * ╚═╝  ╚═╝╚══════╝   ╚═╝   ╚══════╝╚═╝  ╚═╝        

#include <Aster.hpp>
using namespace Aster;

double collect(Graphs::Graph<vec2>* graph, Simulation<vec2>* sim, Body<vec2>* body, Body<vec2>* body2){
    return 1;
}

int main(){
    Composite2d* sim = new Composite2d();
    
    add_composite(sim, 15, 50e7, 50, sim -> get_center(), {0, 0});
    add_composite(sim, 15, 10e6, 50, {sim -> get_width(), sim -> get_height() / 2}, {-.05, .02});
    //add_composite(sim, 15, 10e6, 50, sim -> get_center() * 3/4, {0,0});

    sim 
    -> update_with(LEAPFROG)
    -> get_force_with(PN1)
    -> set_dt(1)
    -> set_scale(12)
    -> load();
    
    render(sim)
    -> show_axis()
    -> show();
}












/*    

    Simulation::presets::add_body(9e11, {WIDTH, 0}, {-5, 3}, sf::Color::Green);
    Simulation::presets::make_ring(HEIGHT/4, HEIGHT/2,{WIDTH/2, HEIGHT/2} ,10000);



    Simulation::presets::add_body(5e11, {WIDTH/8 + WIDTH/9, HEIGHT/2}, {0,  200}, sf::Color::Red);
    Simulation::presets::add_body(5e11, {WIDTH/7 + WIDTH/9, HEIGHT/2}, {0, -200}, sf::Color::Blue);
    Simulation::presets::add_body(9e11, {WIDTH, 0}, {-5, 3}, sf::Color::Green);
*/
