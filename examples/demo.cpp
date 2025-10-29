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

void load_solar_system(Simulation* sim){
    add_body(sim, 1.989e30, {0,0,0}, {0,0, 0}); //sun
    add_body(sim, 3.3011e23, -vec3(0.31 * AU, 0, 0), {0, 38700, 0}); //mercury 
    add_body(sim, 4.8675e24, -vec3(0.71 * AU, 0, 0), {0, 34790, 0}); //venus 
    add_body(sim, 5.9721e24, -vec3(1.00 * AU, 0, 0), {0, 29782, 0}); //earth 
    //add_body(sim, 6.4171e23, -vec3(1.67 * AU, 0, 0), {0, 23130, 0}); //mars
    //add_body(sim, 1.8982e27, -vec3(4.95 * AU, 0, 0), {0, 16060, 0}); //jupiter
    //add_body(sim, 8.6810e25, -vec3(18.3 * AU, 0, 0), {0,  6800, 0}); //uranus
    //add_body(sim, 1.0240e26, -vec3(29.8 * AU, 0, 0), {0,  5430, 0}); //neptune
}
/*
void spinny(Simulation<vec2>* _s, int n, vec2 v, vec2 p){
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

void spinny(Simulation<vec3>* _s, int n, vec3 v, vec3 p){
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
*/
int main(){
    Simulation* _s = new Simulation();
    load_solar_system(_s);
    //add_disk(_s, 1e5, vec3(0), 200, 1, {0,0,0}, 1e-2, 1e30);
    _s -> gravity_solver = GPU_BARNES_HUT;
    _s -> integrator = WH_PLANETARY;
    _s -> integrator_order = 0;
    _s -> force_used = NEWTON;
    
    //_s -> softening = 1;
    _s 
    -> set_scale(200e7) 
    -> set_dt(1e4)
    -> collect_hamiltonian()
    -> collect_distance()
    -> collect_error()
    -> load();

    _s -> integrate(1, true);
    //exit(0);
    
    auto a = Renderer::Renderer3d(_s);
    a.show();
}