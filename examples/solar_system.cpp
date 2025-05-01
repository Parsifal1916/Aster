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
}