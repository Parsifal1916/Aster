// *  █████╗ ███████╗████████╗███████╗██████╗ 
// * ██╔══██╗██╔════╝╚══██╔══╝██╔════╝██╔══██╗
// * ███████║███████╗   ██║   █████╗  ██████╔╝
// * ██╔══██║╚════██║   ██║   ██╔══╝  ██╔══██╗
// * ██║  ██║███████║   ██║   ███████╗██║  ██║
// * ╚═╝  ╚═╝╚══════╝   ╚═╝   ╚══════╝╚═╝  ╚═╝        

#include <Aster.hpp>
#include <string>
#include <cmath>

using namespace Aster;


template <typename F>
void parallel(int cores, size_t num, F func) {
    const size_t n_threads = std::min(num, static_cast<size_t>(cores));
    
    if (num == 0) return;
    
    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    
    const size_t chunk_size = (num + n_threads - 1) / n_threads;
    
    for (size_t i = 0; i < n_threads; ++i) {
        const size_t start = i * chunk_size;
        const size_t end = std::min(start + chunk_size, num);
        
        if (start >= end) break;
        
        threads.emplace_back([start, end, &func]() {
            for (size_t b = start; b < end; ++b) {
                func(b);
            }
        });
    }
    
    for (auto& t : threads) {
        t.join();
    }
}

REAL energy = 0;

template <typename T>
REAL get_total_energy(Simulation<T>* _s){
    REAL retval = 0.0;
    std::mutex retval_mutex;
    const int num_threads = 16;
    int num_bodies = _s->bodies.positions.size();
    
    parallel(_s -> get_cores(), num_bodies, [&retval, &retval_mutex, _s, num_bodies](size_t i) {
        REAL partial_retval = 0.5 * _s -> bodies.get_velocity_of(i).sqr_magn() * _s -> bodies.get_mass_of(i);
        for (int j = i + 1; j < num_bodies; ++j) {
            REAL r = (_s -> bodies.get_position_of(j) - _s -> bodies.get_position_of(i)).magnitude();
            partial_retval -= _s->get_G() * _s -> bodies.get_mass_of(j) * _s -> bodies.get_mass_of(i) / (r);
        }
        
        std::lock_guard<std::mutex> lock(retval_mutex);
        retval += partial_retval;
    });

    return retval;
}

vec2 com = vec2(0);

inline vec2 get_com(Simulation<vec2>* _s){
    vec2 temp = vec2(0);
    std::mutex retval_mutex;
    parallel(_s -> get_cores(), _s -> bodies.positions.size(), [&retval_mutex, &temp, _s](size_t i){
        vec2 partial_retval = _s -> bodies.positions[i] * _s -> bodies.masses[i];

        std::lock_guard<std::mutex> lock(retval_mutex);
        temp += partial_retval;
    });

    return temp /  _s -> get_total_mass();
}

template <typename T>
REAL error_collect(Graphs::Graph<T>* g, Simulation<T>* _s, size_t b){
    com = get_com(_s);
    REAL current_energy = get_total_energy(_s);
    return std::abs(energy - current_energy) / std::abs(energy);
}



REAL get_x(Graphs::Graph<vec2>* g, Simulation<vec2>* _s, size_t b){
    return _s -> bodies.positions[b].x - com.x;
}

REAL get_y(Graphs::Graph<vec2>* g, Simulation<vec2>* _s, size_t b){
    return _s -> bodies.positions[b].y - com.y;
}

int main(){//-> set_theta(.8)
    //Barnes::Barnes_Hut<vec2>* sim = new Barnes::Barnes_Hut<vec2>();

    auto* sim = bake(LIGHT);

    sim 
    -> use_GPU()
    -> set_scale(60e7)
    -> get_force_with(NEWTON)
    -> update_with(SABA10)
    -> set_dt(10e2)
    -> add_graph(error_collect<vec2>)
    -> add_graph(get_x, FOR_EACH)
    -> add_graph(get_y, FOR_EACH)

    ;
 
    double AU = 150e9;
//
    //add_body(sim, 1.989e30, sim -> get_center(), {0,0}); //sole
    //add_body(sim, 3.3011e23, sim -> get_center() - vec2(0.31 * AU, 0), {0, 38700}); //mercury 
    //add_body(sim, 4.8675e24, sim -> get_center() - vec2(0.71 * AU, 0), {0, 34790}); //venus 
    //add_body(sim, 5.9721e24, sim -> get_center() - vec2(1.00 * AU, 0), {0, 29782}); //earth 
    //add_body(sim, 6.4171e23, sim -> get_center() - vec2(1.67 * AU, 0), {0, 23130}); //mars
    //add_body(sim, 1.8982e27, sim -> get_center() - vec2(4.95 * AU, 0), {0, 13060}); //jupiter
    //add_body(sim, 5.6834e26, sim -> get_center() - vec2(9.04 * AU, 0), {0,  9680}); //saturn
    //add_body(sim, 8.6810e25, sim -> get_center() - vec2(18.3 * AU, 0), {0,  6800}); //uranus
    //add_body(sim, 1.0240e26, sim -> get_center() - vec2(29.8 * AU, 0), {0,  5430}); //neptune

    //disk(sim, 10e3, sim -> get_center(), 350, 10, 10e3);
    //cosmic_web(sim, 10e3, 10e2);

    REAL r = 150e8;
    REAL vel = std::sqrt(1.989e30 * sim -> get_G() / r) / 2;

    add_body(sim, 1.989e30, sim -> get_center() - vec2({13*r, 0}), {0,5*vel}); //sole
    add_body(sim, 1.989e30, sim -> get_center() - vec2({12*r, 0}), {0,7*vel}); //sole
    add_body(sim, 1.989e32, sim -> get_center(), {vel/10, vel/2}); //sole add_body(sim, 1.989e29, sim -> get_center() - vec2({15*r, 0}), {vel/8, vel/2});
    
    add_body(sim, 1.989e30, sim -> get_center() + vec2({13*r, 0}), {0,-5*vel}); //sole
    add_body(sim, 1.989e30, sim -> get_center() + vec2({12*r, 0}), {0,-7*vel}); //sole

    sim -> load();
    energy = get_total_energy(sim);    

    sim -> integrate(5e6);//17570

    //for (int i = 0; i < 5e6; ++i){
    //    sim-> step();
    //    std::cout << "step number " << i<< "\n";
    //}
    //render(sim) -> show();
}