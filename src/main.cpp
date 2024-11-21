#include <SFML/Graphics.hpp>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <chrono>
#include <iomanip>
#include <thread>
#include <ctime>

double time_passed = 0;

#define spin .9// v -> 1/spin
#define stV 0

#include "physics/vectors.h"
#include "physics/body.h"
#include "physics/physics.h"


#include "sim-setup/sim_helper.h"
#include "sim-setup/presets.h"
#include "sim-setup/simulation.h"

#include "graphics-api/canvas.h"



int main() {
    //std::cout << sizeof(Body);
    //exit(0);
    rng.seed(static_cast<unsigned int>(std::time(0)));

    Simulation::init();
    Simulation::populate();

    using Clock = std::chrono::high_resolution_clock;
    using Duration = std::chrono::duration<double>;
    
    auto p_start = Clock::now();

    sf::Texture texture;
    sf::Sprite sprite;
    auto prev_time = p_start;
    double elapsed = 0;
    setup_canvas();
    
    while (time_passed <= 50000) {
        auto start = Clock::now();

        Simulation::step();
        draw_bodies();
        //draw_circle();

        auto end = Clock::now();
        // Aggiorna la texture e il sprite con l'immagine
        image.saveToFile("file" + std::to_string((int)(time_passed)) + ".png");
        double elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(Clock::now() - p_start).count();
        double estimated_total_seconds = elapsed_seconds * (50000.0 / (time_passed + 1));
        double remaining_seconds = estimated_total_seconds - elapsed_seconds;

        texture.loadFromImage(image);
        sprite.setTexture(texture);
        Duration cl = Clock::now() - p_start;
        std::cout << "[" << cl.count() << "\t] finished and saved frame number " << time_passed <<" of simulation with scale " << Simulation::current_a << " (took " << (end-  p_start).count() << "ms) \n";
        
        time_passed++;
        image = blank;
    }

    return 0;
}
