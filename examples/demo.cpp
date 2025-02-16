// *  █████╗ ███████╗████████╗███████╗██████╗ 
// * ██╔══██╗██╔════╝╚══██╔══╝██╔════╝██╔══██╗
// * ███████║███████╗   ██║   █████╗  ██████╔╝
// * ██╔══██║╚════██║   ██║   ██╔══╝  ██╔══██╗
// * ██║  ██║███████║   ██║   ███████╗██║  ██║
// * ╚═╝  ╚═╝╚══════╝   ╚═╝   ╚══════╝╚═╝  ╚═╝        

#include <Aster.hpp>
using namespace Aster;

#include <vector>
#include <thread>
#include <mutex>
#include <iostream>
#include <cmath>



double collector(Graphs::Graph<vec2>* g, Simulation<vec2>* _s, Body<vec2>* b) {
    double retval = 0.0;
    const int num_threads = 16;
    int num_bodies = _s->bodies.size();
    std::mutex retval_mutex;
    
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
    auto* sim = bake(LIGHT); //Composite2d* sim = new Composite2d();
    
    sim 
    -> update_with(SYMPLECTIC4)
    -> get_force_with(NEWTON)
    -> set_dt(1e15)
    -> set_scale(1.49e11)
    -> add_graph(collector, ONCE)
    -> load();
    
    add_body(sim, 2e39, sim -> get_center(), {0,0});
    add_body(sim, 5.98e24, {sim -> get_width(), sim -> get_height()/2}, {0, 3e-4});

    //add_composite(sim, 15, 10e6, 50, {0, sim -> get_height() / 2}, {.01, -.002}, 1);
    //add_composite(sim, 15, 10e6, 50, {sim -> get_width(), sim -> get_height() / 2}, {-.01, .002}, 1);
    //add_composite(sim, 15, 10e6, 50, sim -> get_center() * 3/4, {0,0});
    //add_disk(sim, 1e4, sim -> get_center(), 300, 10);
    
    render(sim) -> show();
}












/*    
int main(){
    Composite2d* sim = new Composite2d();
    
    add_composite(sim, 15, 10e6, 50, {0, sim -> get_height() / 2}, {.01, -.002}, 1);
    add_composite(sim, 15, 10e6, 50, {sim -> get_width(), sim -> get_height() / 2}, {-.01, .002}, 1);
    //add_composite(sim, 15, 10e6, 50, sim -> get_center() * 3/4, {0,0});

    sim 
    -> update_with(LEAPFROG)
    -> get_force_with(PN1)
    -> set_dt(10)
    -> set_scale(2)
    -> load();
    
    render(sim)
    -> show_axis()
    -> show();
}
*/
