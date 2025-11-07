// *  █████╗ ███████╗████████╗███████╗██████╗ 
// * ██╔══██╗██╔════╝╚══██╔══╝██╔════╝██╔══██╗
// * ███████║███████╗   ██║   █████╗  ██████╔╝
// * ██╔══██║╚════██║   ██║   ██╔══╝  ██╔══██╗
// * ██║  ██║███████║   ██║   ███████╗██║  ██║
// * ╚═╝  ╚═╝╚══════╝   ╚═╝   ╚══════╝╚═╝  ╚═╝        

#include <Aster.hpp>
#include <string>
using namespace Aster;
#include <math.h>
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

void gravitating_sphere(Simulation* _s, size_t nums, vec3 center, REAL inner, REAL outer, REAL avr_mass, REAL center_mass, vec3 v = {0,0,0}){
    if (warn_if(!_s, "the given simulation to make the gravitating_sphere is a nullptr"))
        return;

    if (warn_if(!nums, "cannot generate a gravitating_sphere with fewer than zero bodies"))
        return;

    if (warn_if(avr_mass <= 0, "gravitating_sphere cannot have less than zero mass"))
        return;

    if (warn_if(inner <= 0, "gravitating_sphere cannot have negative radius"))
        return;


    // sets up the cluster
    Cluster cluster;
    cluster.number = nums;
    cluster.name = "gravitating_sphere";
    add_body(_s, 
        center_mass, 
        center,
        v
    );

    // sets up the builder lambda
    cluster.builder = [inner, outer, v, avr_mass, center, _s, center_mass ](Cluster cl3d, size_t _) {
        REAL r = rng_val(inner, outer);
        REAL theta = rng_val(0, M_PI);
        REAL phi = rng_val(0, 2*M_PI);

        vec3 pos = {
            r * sin(theta) * cos(phi),
            r * sin(theta) * sin(phi),
            r * cos(theta)
        };

        vec3 r_hat = pos.normalize();

        vec3 arbitrary = (fabs(r_hat.x) < 0.9) ? vec3{1,0,0} : vec3{0,1,0};
        vec3 t1 = r_hat.cross(arbitrary).normalize();
        vec3 t2 = r_hat.cross(t1);

        REAL angle = rng_val(0, 2*M_PI);
        vec3 tangent_dir = cos(angle) * t1 + sin(angle) * t2;

        REAL magn_vel = std::sqrt(_s->get_G() * center_mass / r);
        vec3 vel = tangent_dir * magn_vel;

        add_body(_s,
            rng_percent() * avr_mass,
            pos + center,
            vel + v
        );
    };

    // adds the cluster to the queue
    _s -> loading_queue.add_cluster(cluster);
}

int main() {
    Simulation* sim = new Simulation();
    add_disk(sim, 1e2, vec3(0), 200, 1, {0,0,0}, 1e-2, 1e30);
    
    sim->gravity_solver = BARNES_HUT;
    sim->integrator = WH_PLANETARY;
    sim->integrator_order = 0;
    sim->force_used = NEWTON;
    sim -> theta = .44;
    sim->softening = 1e-16;
    sim->set_scale(200e7)->set_dt(1E4)->load();
    auto report = sim->integrate(3, true).time_per_step;
}
