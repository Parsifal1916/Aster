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
    add_body(sim, 3.3011e23, sim -> get_center() - vec3(0.31 * AU, 0, 0), {0, 38700, 0}); //mercury 
    add_body(sim, 1.989e30, sim -> get_center(), {0,0, 0}); //sun
    add_body(sim, 4.8675e24, sim -> get_center() - vec3(0.71 * AU, 0, 0), {0, 34790, 0}); //venus 
    add_body(sim, 5.9721e24, sim -> get_center() - vec3(1.00 * AU, 0, 0), {0, 29782, 0}); //earth 
    add_body(sim, 6.4171e23, sim -> get_center() - vec3(1.67 * AU, 0, 0), {0, 23130, 0}); //mars
    add_body(sim, 1.8982e27, sim -> get_center() - vec3(4.95 * AU, 0, 0), {0, 16060, 0}); //jupiter
    //add_body(sim, 1.0000e30, sim -> get_center() + v3c2(9.04 * AU,, 0 0), {0, -96, 080}); //saturn
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

template <typename type>
void test_sim(int n){
    Barnes::BH_hyper<type>* sim = new Barnes::BH_hyper<type>;
    add_disk(sim, n, sim -> get_center(), 200, 10, 10e13);
    
    sim 
    -> set_theta(1)
    -> set_dt(1)
    -> set_scale(200e7)
    -> get_force_with(NEWTON)
    -> update_with(SABA1)
    -> load()
    ;
    

    sim -> load();

    sim -> integrate(10);
}

void spinny(Barnes::BH_hyper<vec2>* _s, int n, vec2 v, vec2 p){
    double M = 10e34;
    add_body(_s, M, p, v);

    Cluster<vec2> cluster;
    cluster.number = n;
    cluster.name = "Disk";

    // creates the builder lambda
    cluster.builder = [v, M, p, _s ](Cluster<vec2> cl2d, size_t _) {
        vec2 pos = rng_point_in_circle(100, 30)* _s -> get_scale(); // gets a random point inside the disk
        //pos.y /= 2.5;
        // generates the radius from the position
        REAL radius = pos.magnitude();

        // velocty on that point
        REAL magn_vel = std::sqrt(_s -> get_G() * M/radius);
        vec2 vel = vec2(pos.y/radius * magn_vel, -pos.x/radius * magn_vel);
       

        // assembles the body
        add_body(_s,
            rng_percent() * 2e20, 
            pos + p,
            vel + v
        );
    };

    // adds the cluster to the queue
    _s -> loading_queue.add_cluster(cluster);
}

void spinny(Barnes::BH_hyper<vec3>* _s, int n, vec3 v, vec3 p){
    double M = 10e34;
    add_body(_s, M, p, v);

    Cluster<vec3> cluster;
    cluster.number = n;
    cluster.name = "Disk";

    // creates the builder lambda
    cluster.builder = [v, M, p, _s ](Cluster<vec3> cl2d, size_t _) {
        vec3 pos = rng_point_in_cylinder(100, 30)* _s -> get_scale(); // gets a random point inside the disk
        //pos.y /= 2.5;
        // generates the radius from the position
        REAL radius = pos.magnitude();

        // velocty on that point
        REAL magn_vel = std::sqrt(_s -> get_G() * M/radius);
        vec3 vel = vec3(pos.y/radius * magn_vel, -pos.x/radius * magn_vel, 0);
       

        // assembles the body
        add_body(_s,
            rng_percent() * 2e20, 
            pos + p,
            vel + v
        );
    };

    // adds the cluster to the queue
    _s -> loading_queue.add_cluster(cluster);
}

int main(){ 
    using namespace Aster::Graphs;
    size_t 
        bodies = 5000,
        steps = 1000
    ;  
    Barnes::BH_hyper<vec2>* sim = new Barnes::BH_hyper<vec2>;
    //load_solar_system(sim);
    //add_disk(sim, bodies, sim -> get_center(), 200, 10, {0,0,0}, 1, {0,0,0});
    //spinny(sim, bodies, {0,0,0}, {0, 0,0});
    //cosmic_web(sim, bodies, bodies);
    //add_body(sim, 50e32, {sim -> get_width()/3, 0}, {0, 5e6});
    spinny(sim, bodies, {0,0}, {0, 0});
    //cosmic_web(sim, 2e5, 2e27);


    //std::vector<vec2> pos, vel;
    //std::vector<double> mass;

    sim 
    -> always_read_positions()
    -> set_theta(1)
    -> set_dt(1e4)
    -> set_scale(200e7)
    -> get_force_with(NEWTON)
    -> update_with(SABA1)
    -> load()
    ;

    //pos = sim -> bodies.positions;
    //vel = sim -> bodies.velocities;
    //mass = sim -> bodies.masses;
    render(sim) -> show();exit(0);

    double before =  sim -> get_total_energy();
    sim -> integrate(steps);
    sim -> read_positions();
    

    double after = sim -> get_total_energy();
    std::cout << "sim1 energy before = " << before << "\n";
    std::cout << "sim1 energy after  = " << after << "\n";
    std::cout << "sim1 ratio  = " << std::abs(after  - before)/ before * 10e5 << "\n";
    Barnes::Barnes_Hut<vec3>* sim2 = new Barnes::Barnes_Hut<vec3>;
    load_solar_system(sim2);
    //add_disk(sim2, bodies, sim -> get_center(), 800, 10, 10e27);

    //for (int i = 0; i < mass.size(); ++i){
    //    add_body(sim2, mass[i], pos[i], vel[i]);
    //}

    sim2
    -> set_theta(1)
    -> set_dt(1e4)
    -> set_scale(200e7)
    -> get_force_with(NEWTON)
    -> update_with(SABA1)
    -> load()
    ;
    
    before = get_total_energy(sim2);
    
    sim2 -> integrate(steps);

    std::cout << "sim2 energy before = " << before << "\n";
    after = get_total_energy(sim2);
    std::cout << "sim2 energy after  = " << after << "\n";
    std::cout << "sim2 ratio  = " << std::abs(after - before)/after * 10e5;
}













/*    

    Simulation::presets::add_body(9e11, {WIDTH, 0}, {-5, 3}, sf::Color::Green);
    Simulation::presets::make_ring(HEIGHT/4, HEIGHT/2,{WIDTH/2, HEIGHT/2} ,10000);



    Simulation::presets::add_body(5e11, {WIDTH/8 + WIDTH/9, HEIGHT/2}, {0,  200}, sf::Color::Red);
    Simulation::presets::add_body(5e11, {WIDTH/7 + WIDTH/9, HEIGHT/2}, {0, -200}, sf::Color::Blue);
    Simulation::presets::add_body(9e11, {WIDTH, 0}, {-5, 3}, sf::Color::Green);
*/
