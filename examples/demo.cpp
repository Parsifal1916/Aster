// *  █████╗ ███████╗████████╗███████╗██████╗ 
// * ██╔══██╗██╔════╝╚══██╔══╝██╔════╝██╔══██╗
// * ███████║███████╗   ██║   █████╗  ██████╔╝
// * ██╔══██║╚════██║   ██║   ██╔══╝  ██╔══██╗
// * ██║  ██║███████║   ██║   ███████╗██║  ██║
// * ╚═╝  ╚═╝╚══════╝   ╚═╝   ╚══════╝╚═╝  ╚═╝        

#include <Aster.hpp>
#include <string>
using namespace Aster;

#define AU 150e9

void load_solar_system(Simulation<vec3>* sim){
    add_body(sim, 1.989e30, sim -> get_center(), {0,0,0}); //sun
    add_body(sim, 3.3011e23, sim -> get_center() - vec3(0.31 * AU, 0, 0), {0, 38700, 0}); //mercury 
    add_body(sim, 4.8675e24, sim -> get_center() - vec3(0.71 * AU, 0, 0), {0, 34790, 0}); //venus 
    add_body(sim, 5.9721e24, sim -> get_center() - vec3(1.00 * AU, 0, 0), {0, 29782, 0}); //earth 
    add_body(sim, 6.4171e23, sim -> get_center() - vec3(1.67 * AU, 0, 0), {0, 23130, 0}); //mars
    add_body(sim, 1.8982e27, sim -> get_center() - vec3(4.95 * AU, 0, 0), {0, 16060, 0}); //jupiter
    //add_body(sim, 1.0000e30, sim -> get_center() + vec3(9.04 * AU, 0, 0), {0, -9680, 0}); //saturn
    add_body(sim, 8.6810e25, sim -> get_center() - vec3(18.3 * AU, 0, 0), {0,  6800, 0}); //uranus
    add_body(sim, 1.0240e26, sim -> get_center() - vec3(29.8 * AU, 0, 0), {0,  5430, 0}); //neptune
}//

void load_solar_system(Simulation<vec2>* sim){
    add_body(sim, 3.3011e23, sim -> get_center() - vec2(0.31 * AU, 0), {0, 38700}); //mercury 
    add_body(sim, 1.989e30, sim -> get_center(), {0,0}); //sun
    add_body(sim, 4.8675e24, sim -> get_center() - vec2(0.71 * AU, 0), {0, 34790}); //venus 
    add_body(sim, 5.9721e24, sim -> get_center() - vec2(1.00 * AU, 0), {0, 29782}); //earth 
    add_body(sim, 6.4171e23, sim -> get_center() - vec2(1.67 * AU, 0), {0, 23130}); //mars
    add_body(sim, 1.8982e27, sim -> get_center() - vec2(4.95 * AU, 0), {0, 16060}); //jupiter
    //add_body(sim, 1.0000e30, sim -> get_center() + vec2(9.04 * AU, 0), {0, -9680}); //saturn
    add_body(sim, 8.6810e25, sim -> get_center() - vec2(18.3 * AU, 0), {0,  6800}); //uranus
    add_body(sim, 1.0240e26, sim -> get_center() - vec2(29.8 * AU, 0), {0,  5430}); //neptune
}//

template <typename T>
void update_sim(Simulation<T>* _s){
    for (int i = 0; i < _s -> bodies.positions.size(); i++){
        //std::cout << _s -> bodies.get_acc_of(i).sqr_magn() << ", " << _s -> bodies.get_position_of(i).sqr_magn() << "\n";
        _s -> bodies.get_velocity_of(i) += _s -> bodies.get_acc_of(i) * _s -> get_dt();
        _s -> bodies.get_position_of(i) += _s -> bodies.get_velocity_of(i) * _s -> get_dt();
    }
}

void add_bodies(Simulation<vec3>* _s){
    int outer = 100, inner = 10;
    double avr_mass = 10e16;
    int N = 1e4;
    for (int i = 0; i < N; ++i){
        vec3 pos = rng_point_in_cylinder(outer, inner, 0)* _s -> get_scale(); // get*s a random point inside the disk

        // generates the radius from the position
        double radius = pos.sqr_magn();

        // velocty on that point
        double magn_vel =std::sqrt(_s -> get_G() * avr_mass * N / radius);
        vec3 vel = vec3(-pos.y/radius * magn_vel, pos.x/radius * magn_vel, 0);

        // assembles the body
        add_body(_s,
            avr_mass, 
            pos + _s ->get_center(),
            vel
        );
    }
}

//template <typename type>
//void test_sim(float n){
//    Barnes::BHG<type>* sim = new Barnes::BHG <type>;
//
//    sim 
//    -> set_theta(0.8)
//    -> use_GPU()
//    -> set_scale(1000)
//    -> set_dt(1)
//    -> update_with(update_sim<type>);
//    ;
//    add_disk(sim, n, sim -> get_center(), 200, 10, 10e13);
//
//    sim -> load();
//
//    sim -> integrate(10);
//    delete sim;
//}

int main(){ 
    //for (const auto& i : {100, 150, 500, 1000, 1500, 10000, 50000, 100000, 500000, 1000000}){
    //    test_sim<vec2>(i);
    //}
    using type = vec2;

    auto* sim = new Barnes::BHG <type>;

    sim 
    -> set_theta(.9)
    //-> use_GPU()
    -> update_with(SABA6)
    -> set_scale(150e6)
    -> set_dt(10e3)
    -> update_with(update_sim<type>);
    ;
    //load_solar_system(sim);
    add_disk(sim, 1e6, sim -> get_center(), 200, 10, 10e16);

    sim -> load();
    //sim -> step();
    //render(sim) -> show();
    sim -> integrate(100);
}













/*    

    Simulation::presets::add_body(9e11, {WIDTH, 0}, {-5, 3}, sf::Color::Green);
    Simulation::presets::make_ring(HEIGHT/4, HEIGHT/2,{WIDTH/2, HEIGHT/2} ,10000);



    Simulation::presets::add_body(5e11, {WIDTH/8 + WIDTH/9, HEIGHT/2}, {0,  200}, sf::Color::Red);
    Simulation::presets::add_body(5e11, {WIDTH/7 + WIDTH/9, HEIGHT/2}, {0, -200}, sf::Color::Blue);
    Simulation::presets::add_body(9e11, {WIDTH, 0}, {-5, 3}, sf::Color::Green);
*/
