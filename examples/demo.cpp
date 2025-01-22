// *  █████╗ ███████╗████████╗███████╗██████╗ 
// * ██╔══██╗██╔════╝╚══██╔══╝██╔════╝██╔══██╗
// * ███████║███████╗   ██║   █████╗  ██████╔╝
// * ██╔══██║╚════██║   ██║   ██╔══╝  ██╔══██╗
// * ██║  ██║███████║   ██║   ███████╗██║  ██║
// * ╚═╝  ╚═╝╚══════╝   ╚═╝   ╚══════╝╚═╝  ╚═╝        

#include <Aster.hpp>
using namespace Aster;

int main(){
    auto* sim = bake(BARNES_HUT);
    cosmic_web(sim, 1e4, 10e10);
    //add_disk3d(sim, 1e4, sim -> get_center(), 60e1, .3, {}, 10e10);
    //add_disk3d(sim, 1e4, sim -> get_corner(0), 60e1, .3, {90, 45, 0}, 10e10);
    //add_disk(sim, 1e3, sim -> get_center(), 100, 10);
    //rng_sphere(sim, 10e3, {meta.WIDTH, meta.HEIGHT, meta.depth}, 10e2);
    
    sim 
    -> set_dt(0)
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
