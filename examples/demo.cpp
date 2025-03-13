// *  █████╗ ███████╗████████╗███████╗██████╗ 
// * ██╔══██╗██╔════╝╚══██╔══╝██╔════╝██╔══██╗
// * ███████║███████╗   ██║   █████╗  ██████╔╝
// * ██╔══██║╚════██║   ██║   ██╔══╝  ██╔══██╗
// * ██║  ██║███████║   ██║   ███████╗██║  ██║
// * ╚═╝  ╚═╝╚══════╝   ╚═╝   ╚══════╝╚═╝  ╚═╝        

#include <Aster.hpp>
using namespace Aster;

double collect_x(Graphs::Graph<vec2>* g, Simulation<vec2>* _s, Body<vec2>* b){
    vec2 bar = {0,0};
    double mass = 0;


    for (const auto& body : _s -> bodies){
        mass += body.mass;
        bar += body.position*body.mass;
    }

    bar = bar / mass;
    return bar.x - b -> position.x;
}

double collect_y(Graphs::Graph<vec2>* g, Simulation<vec2>* _s, Body<vec2>* b){
    vec2 bar = {0,0};
    double mass = 0;


    for (const auto& body : _s -> bodies){
        mass += body.mass;
        bar += body.position*body.mass;
    }

    bar = bar / mass;
    return bar.y - b -> position.y;
}

int main(){
    auto* sim = bake(LIGHT);


    sim 
    -> set_scale(10e10)
    -> get_force_with(NEWTON)
    -> set_dt(2e8)
    -> add_graph(collect_x, FOR_EACH)
    -> add_graph(collect_y, FOR_EACH)
    -> load();
 
    add_body(sim, 2e30, sim -> get_center(), {0,0});
    add_body(sim, 2e14, {sim -> get_width() *6/7,  sim -> get_height() - 10e10*768 /2 },  vec2({0, 600}));
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
