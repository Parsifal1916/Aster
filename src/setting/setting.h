#pragma once

#include "../sim-setup/simulation.h"
#include "../sim-setup/sim_helper.h"
#include <cassert>

namespace Simulation{
    using Clock = std::chrono::high_resolution_clock;
    using Duration = std::chrono::duration<double>;

    void set_numof_objs(unsigned int n_){
        obj = n_;
        if (obj < NUM_THREADS)
            NUM_THREADS = obj;
    }

    /*
    * set screen's height, width and bg color (optional)
    * @param w_: screen's height
    * @param h_: screen's width
    * @param c_: screen's bg color
    */
    void set_screen_size(unsigned int w_, unsigned int h_, sf::Color c_ = sf::Color::Black){
        HEIGHT = h_;
        WIDTH = w_;
        bg = c_;
    }

    /*
    * sets the maximum number of threads for the given simulation
    * if the given number is higher than the number of objs
    * it will be topped at the number of objs
    */
    void set_max_threads(unsigned int t_){
        NUM_THREADS = t_;
        if (t_ == 0)
            NUM_THREADS = 1;
        if (obj < NUM_THREADS)
            NUM_THREADS = obj;
    }

    void set_dt(float dt_){
        assert(dt_ > 0);
        dt = dt_;
    }

    void set_omega_m(float om_){
        assert(om_ > 0);
        omega_m = om_;
    }

    void set_omega_l(float ol_){
        assert(ol_ > 0);
        omega_l = ol_;
    }

    void set_hubble(float h_){
        assert(h_ > 0);
        H_0 = h_;
    }

    void set_vacuum_density(float d_){
        assert(d_ > 0);
        vacuum_density = d_;
        root_omega = sqrt(d_);
    }

    void set_max_frames(unsigned int f_){
        max_frames = f_;
    }

    void set_sim_type(short type){
        assert(type < 4);
        sim_type = static_cast<simulation_types>(type);
    }

    // !! deprecated
    void save_stebystep(){
        auto p_start = Clock::now();

        sf::Texture texture;
        sf::Sprite sprite;
        auto prev_time = p_start;
        double elapsed = 0;
        Renderer::setup_canvas();

        while (time_passed <= 50000) {
            auto start = Clock::now();

            Simulation::step();
            Renderer::render();

            auto end = Clock::now();
           
            double elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(Clock::now() - p_start).count();
            double estimated_total_seconds = elapsed_seconds * (50000.0 / (time_passed + 1));
            double remaining_seconds = estimated_total_seconds - elapsed_seconds;

            texture.loadFromImage(Renderer::image);
            sprite.setTexture(texture);
            Duration cl = Clock::now() - p_start;
            std::cout << "[" << cl.count() << "\t] finished and saved frame number " << time_passed <<" of simulation with scale " << Simulation::current_a << " (took " << (end-  p_start).count() << "ms) \n";
            
            time_passed++;
            Renderer::reset();
        }
    }

    /*
    * spawn a window and starts the simulation
    */
    void do_window(){
        sf::Texture texture;
        sf::Sprite sprite;

        // sets up a window and a texture renderer
        Renderer::setup_canvas();
        Renderer::window.create(sf::VideoMode(WIDTH, HEIGHT), "Simulation");
        
        sf::RenderTexture render_txture;
        render_txture.create(WIDTH, HEIGHT);
        
        while (Renderer::window.isOpen()){
            sf::Event event;
            // handles the window closing
            while (Renderer::window.pollEvent(event)){
                if (event.type == sf::Event::Closed) Renderer::window.close();
            }    

            // take a timestep and renders it
            step();
            Renderer::render();
            Renderer::window.clear();

            // overlays the lagrangian text object
            render_txture.clear();
            render_txture.draw(sprite);
            render_txture.draw(Renderer::lagr_text);
            render_txture.display();

            texture.loadFromImage(Renderer::image);
            sprite.setTexture(texture);
            Renderer::window.draw(sprite);

            if (save)
                render_txture.getTexture().copyToImage().saveToFile("file" + std::to_string((int)(time_passed)) + ".png");

            Renderer::window.draw(Renderer::lagr_text);
            Renderer::window.display();

            if (!leave_traces)
                Renderer::reset();

            time_passed++;
        }
    }

    /*
    * stats the simulation
    */
    void fire_up(){
        step = general_step;

        obj = bodies.size();
        rng.seed(static_cast<unsigned int>(std::time(0)));
        set_max_threads(NUM_THREADS);
    
        do_window();
    }

    namespace Sauce{
        void set_divergence(float d_){
            assert(d_ > 0);
            divergence = d_;
        }

        void set_takeover(float t_){
            assert(t_ > 0);
            takeover = t_;
        }
    }
}