// *  █████╗ ███████╗████████╗███████╗██████╗ 
// * ██╔══██╗██╔════╝╚══██╔══╝██╔════╝██╔══██╗
// * ███████║███████╗   ██║   █████╗  ██████╔╝
// * ██╔══██║╚════██║   ██║   ██╔══╝  ██╔══██╗
// * ██║  ██║███████║   ██║   ███████╗██║  ██║
// * ╚═╝  ╚═╝╚══════╝   ╚═╝   ╚══════╝╚═╝  ╚═╝        

#include <Aster.hpp>
using namespace Aster;

double collector(Graphs::Graph<vec2>* g, Simulation<vec2>* _s, Body<vec2>* b) {
    double retval = 0.0;
    std::mutex retval_mutex;
    const int num_threads = 16;
    int num_bodies = _s->bodies.size();
    
    auto thread_func = [&retval, &retval_mutex, _s, num_bodies](int start_idx, int end_idx) {
        double partial_retval = 0.0;

        for (int i = start_idx; i < end_idx; ++i) {
            const auto& b1 = _s->bodies[i];
            partial_retval += 0.5 * b1.velocity.sqr_magn() * b1.mass;
            for (int j = i + 1; j < num_bodies; ++j) {
                const auto& b2 = _s->bodies[j];
                double r = (b2.position - b1.position).magnitude();
                partial_retval -= _s->get_G() * b1.mass * b2.mass / r;
            }
        }

        std::lock_guard<std::mutex> lock(retval_mutex);
        retval += partial_retval;
    };

    std::vector<std::thread> threads;
    int chunk_size = num_bodies / num_threads;

    for (int i = 0; i < num_threads; ++i) {
        int start_idx = i * chunk_size;
        int end_idx = (i == num_threads - 1) ? num_bodies : (i + 1) * chunk_size;
        threads.push_back(std::thread(thread_func, start_idx, end_idx));
    }

    for (auto& t : threads) {
        t.join();
    }

    return retval;
}

int main(){
    auto* sim = bake(LIGHT);

    vec2 velocity = vec2({600, -400})*3;

    sim 
    -> set_scale(10e11)
    -> get_force_with(NEWTON)
    -> set_dt(1e8)
    -> load();
 
    add_body(sim, 2e30, {10e10*1366 * 1/4, sim -> get_height() -  10e10*768 /2 }, vec2({0, 800}) + velocity);
    add_body(sim, 2e30, {10e10*1366 *3/4,  sim -> get_height() - 10e10*768 /2 },  vec2({0,-800}) + velocity);

    add_body(sim, 4e31, {sim -> get_width() * 3/4, sim -> get_height() * 1/4}, -velocity + vec2({0, 550}));

    render(sim) -> show();

    
    /*
    auto* sim = bake(BARNES_HUT);
    
    sim 
    -> update_with(EULER)
    -> get_force_with(NEWTON)
    -> set_dt(1e-50)
    -> set_scale(1.49e11 / 4)
    ;
    
    //add_body(sim, 2e30, sim -> get_center(), {0,0});
    //add_body(sim, 5e24, {sim -> get_width() * 3/4, sim -> get_height()/2},{0, 2000});
    //add_body(sim, 5e24, {sim -> get_width(), sim -> get_height()/2},{0, 2790});
    //add_disk3d(sim, 1e4, sim -> get_center(), 60e1, .3, {}, 10e10);
    //add_disk3d(sim, 1e4, sim -> get_corner(0), 60e1, .3, {90, 45, 0}, 10e10);
    add_disk(sim, 1e3, sim -> get_center(), 100, 0);
    sim -> load();
    //rng_sphere(sim, 10e3, {meta.WIDTH, meta.HEIGHT, meta.depth}, 10e2);
    


    render(sim)
    -> show_axis()
    -> show();*/
}













/*    

    Simulation::presets::add_body(9e11, {WIDTH, 0}, {-5, 3}, sf::Color::Green);
    Simulation::presets::make_ring(HEIGHT/4, HEIGHT/2,{WIDTH/2, HEIGHT/2} ,10000);



    Simulation::presets::add_body(5e11, {WIDTH/8 + WIDTH/9, HEIGHT/2}, {0,  200}, sf::Color::Red);
    Simulation::presets::add_body(5e11, {WIDTH/7 + WIDTH/9, HEIGHT/2}, {0, -200}, sf::Color::Blue);
    Simulation::presets::add_body(9e11, {WIDTH, 0}, {-5, 3}, sf::Color::Green);
*/
