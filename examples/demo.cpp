// *  █████╗ ███████╗████████╗███████╗██████╗ 
// * ██╔══██╗██╔════╝╚══██╔══╝██╔════╝██╔══██╗
// * ███████║███████╗   ██║   █████╗  ██████╔╝
// * ██╔══██║╚════██║   ██║   ██╔══╝  ██╔══██╗
// * ██║  ██║███████║   ██║   ███████╗██║  ██║
// * ╚═╝  ╚═╝╚══════╝   ╚═╝   ╚══════╝╚═╝  ╚═╝        

#include <Aster.hpp>
using namespace Aster;

int main(){
    auto meta = sim3d_meta();
    meta.dt = 0;
    meta.type = BARNES_HUT;
    //meta.save = true;

    Simulation3d* sim = bake3d(meta);
    //presets::cosmic_web(sim, 1e4, 10e9);
    ////presets::cosmic_web3d(sim, 10e3, 10e5);
    presets::add_disk3d(sim, 60e1, 10, {meta.WIDTH/2, meta.HEIGHT/2, meta.depth/2}, 1e4, .3, 0, 90, 10e10);
    presets::add_disk3d(sim, 60e1, 10, {0,0,0}, 1e4, .3, 0, 0, 10e10); 
    //presets::rng_sphere(sim, 10e3, {meta.WIDTH, meta.HEIGHT, meta.depth}, 10e2);

    auto n = Renderer::Renderer3d(sim);
    n.run();
}












/*    

    Simulation::presets::add_body(9e11, {WIDTH, 0}, {-5, 3}, sf::Color::Green);
    Simulation::presets::make_ring(HEIGHT/4, HEIGHT/2,{WIDTH/2, HEIGHT/2} ,10000);



    Simulation::presets::add_body(5e11, {WIDTH/8 + WIDTH/9, HEIGHT/2}, {0,  200}, sf::Color::Red);
    Simulation::presets::add_body(5e11, {WIDTH/7 + WIDTH/9, HEIGHT/2}, {0, -200}, sf::Color::Blue);
    Simulation::presets::add_body(9e11, {WIDTH, 0}, {-5, 3}, sf::Color::Green);
*/
