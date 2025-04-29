// *  █████╗ ███████╗████████╗███████╗██████╗ 
// * ██╔══██╗██╔════╝╚══██╔══╝██╔════╝██╔══██╗
// * ███████║███████╗   ██║   █████╗  ██████╔╝
// * ██╔══██║╚════██║   ██║   ██╔══╝  ██╔══██╗
// * ██║  ██║███████║   ██║   ███████╗██║  ██║
// * ╚═╝  ╚═╝╚══════╝   ╚═╝   ╚══════╝╚═╝  ╚═╝        

#include <Aster.hpp>
#include <string>
using namespace Aster;


int main(){
    auto* sim = bake(LIGHT);

    sim 
    -> set_scale(150e7)
    -> get_force_with(NEWTON)
    -> update_with(SABA5)
    -> set_dt(10e3)
    -> collect_distance()
    -> collect_error()
    -> collect_hamiltonian()
    ;
 
    double AU = 150e9;

    add_body(sim, {1.989e30, sim -> get_center(), {0,0}}); //sole
    add_body(sim, {3.3011e23, sim -> get_center() - vec2(0.31 * AU, 0), {0, 38700}}); //mercury 
    add_body(sim, {4.8675e24, sim -> get_center() - vec2(0.71 * AU, 0), {0, 34790}}); //venus 
    add_body(sim, {5.9721e24, sim -> get_center() - vec2(1.00 * AU, 0), {0, 29782}}); //earth 
    add_body(sim, {6.4171e23, sim -> get_center() - vec2(1.67 * AU, 0), {0, 23130}}); //mars
    add_body(sim, {1.8982e27, sim -> get_center() - vec2(4.95 * AU, 0), {0, 13060}}); //jupiter
    add_body(sim, {5.6834e26, sim -> get_center() - vec2(9.04 * AU, 0), {0,  9680}}); //saturn
    add_body(sim, {8.6810e25, sim -> get_center() - vec2(18.3 * AU, 0), {0,  6800}}); //uranus
    add_body(sim, {1.0240e26, sim -> get_center() - vec2(29.8 * AU, 0), {0,  5430}}); //neptune

    sim -> load();

    render(sim) -> show();

    /*
    auto* sim = bake3d(BARNES_HUT);
    
    sim 
    -> update_with(EULER)
    -> get_force_with(NEWTON)
    -> set_dt(1e-50)
    -> set_scale(1)
    ;
    
    //add_body(sim, 2e30, sim -> get_center(), {0,0});
    //add_body(sim, 5e24, {sim -> get_width() * 3/4, sim -> get_height()/2},{0, 2000});
    //add_body(sim, 5e24, {sim -> get_width(), sim -> get_height()/2},{0, 2790});
    //add_disk3d(sim, 1e4, sim -> get_center(), 60e1, .3, {}, 10e10);
    add_disk(sim, 1e4, sim -> get_center(), 60e1, .3, {90, 45, 0}, 10e10);
    //add_disk(sim, 1e3, sim -> get_center(), 100, 0);
    //sim -> load();
    //rng_sphere(sim, 10e3, {meta.WIDTH, meta.HEIGHT, meta.depth}, 10e2);
    


    render(sim)
    -> show_axis()
    -> show();
*/}













/*    

    Simulation::presets::add_body(9e11, {WIDTH, 0}, {-5, 3}, sf::Color::Green);
    Simulation::presets::make_ring(HEIGHT/4, HEIGHT/2,{WIDTH/2, HEIGHT/2} ,10000);



    Simulation::presets::add_body(5e11, {WIDTH/8 + WIDTH/9, HEIGHT/2}, {0,  200}, sf::Color::Red);
    Simulation::presets::add_body(5e11, {WIDTH/7 + WIDTH/9, HEIGHT/2}, {0, -200}, sf::Color::Blue);
    Simulation::presets::add_body(9e11, {WIDTH, 0}, {-5, 3}, sf::Color::Green);
*/
