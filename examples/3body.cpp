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
    return b -> position.x -bar.x;
}

double collect_y(Graphs::Graph<vec2>* g, Simulation<vec2>* _s, Body<vec2>* b){
    vec2 bar = {0,0};
    double mass = 0;


    for (const auto& body : _s -> bodies){
        mass += body.mass;
        bar += body.position*body.mass;
    }

    bar = bar / mass;
    return b -> position.y - bar.y;
}

int main(){
    auto* sim = bake(LIGHT);

    std::cout << "started";

    vec2 velocity = vec2({60, -40})*3;

    sim 
    -> set_scale(10e11)
    -> get_force_with(PN1)
    -> set_dt(2e8)
    -> add_graph(collect_x, FOR_EACH)
    -> add_graph(collect_y, FOR_EACH)
    -> load();
 
    add_body(sim, 2e30, {10e10*1366 * 1/4, sim -> get_height() -  10e10*768 /2 }, vec2({0, 800}) + velocity);
    add_body(sim, 2e30, {10e10*1366 *3/4,  sim -> get_height() - 10e10*768 /2 },  vec2({0,-800}) + velocity);

    add_body(sim, 5e30, {sim -> get_width() * 3/4, sim -> get_height() * 1/4}, -velocity + vec2({0, 170}));

    //render(sim) -> show();

    for (int i = 0; i < 16e3; i++)
        sim -> step();
}

