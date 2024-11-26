#pragma once
#include <cassert>
#include <SFML/Graphics.hpp>

namespace Simulation{
namespace Renderer{
    using render_func = void(*)();

    /*
    * this file contains all the functions
    * relative to image inits, drawing and
    * graphic related functions 
    */

    sf::Image blank;
    sf::Image image;
    sf::Font font;
    sf::Text lagr_text;
    sf::RenderWindow window;
    
    /*
    * initializes the rendering module
    */
    void init(){
        font.loadFromFile("/usr/share/fonts/TTF/Hack-Regular.ttf");
        lagr_text.setFont(font);
        lagr_text.setString("Langragian: ");
        lagr_text.setCharacterSize(16);
        lagr_text.setFillColor(sf::Color::White);
    }
    

    /*
    * cleans the image and gets called every
    * time the simulator want to clear the
    * screen
    */
    void reset(){
        image = blank;
    }

    /*
    * initializes an image reference 
    * as a blank image with a Simulation::bg 
    * background by setting all pixels to that color
    * @param _i: reference to the image
    */
    void init_image(sf::Image& _i){
        _i.create(WIDTH, HEIGHT, bg);
        for (unsigned int x = 0; x < WIDTH; ++x) {
            for (unsigned int y = 0; y < HEIGHT; ++y) {
                _i.setPixel(x, y, bg);
            }
        }
    }


    /*
    * draws a circle of radius $radius and center = (centerX, centerY)
    * note: it draws it on the "image" object
    * @param radius: radius of the circle
    * @param centerX: x component
    * @param centerY: Y component
    */
    void draw_circle(int radius, vec2 center, sf::Color color = sf::Color::White){
        for (int y = center.y - radius; y <= center.y + radius; ++y) {
            for (int x = center.x - radius; x <= center.x + radius; ++x) {
                int dx = x - center.x;
                int dy = y - center.y;
                if (dx * dx + dy * dy <= radius * radius) {
                    if (x >= 0 && x < WIDTH && y >= 0 && y < HEIGHT)
                        image.setPixel(x, y, color);
                }
            }
        }
    }

    /*
    * sets up the "image" object by
    * setting every pixel to the Simulation::bg
    * color
    */
    void setup_canvas(){
        blank.create(WIDTH, HEIGHT, bg);
        for (unsigned int x = 0; x < WIDTH; ++x) {
                for (unsigned int y = 0; y < HEIGHT; ++y) {
                        blank.setPixel(x, y, bg);
                }
        }

        image.create(WIDTH, HEIGHT, bg);
        image = blank;
    }


    void clear_graph(){
        for (int i = 0; i < WIDTH; i++){
        for (int j = graph_height; j < HEIGHT; j++){
            image.setPixel(i, j, bg);
        }
        }
    }

    std::string to_notation(double n){
        std::ostringstream scien;
        scien << std::scientific << std::setprecision(3) << n;
        return scien.str();
    }

    /*
    * draws the langragian graph on the bottom of the screen
    */
    void draw_graph(){
        // sets up and renders the text
        lagr_text.setString("Lagrangian: "+ to_notation(lagrangians.back()) + "J");
        lagr_text.setPosition(10, graph_height + 10);

        window.draw(lagr_text);

        // draws an horizontal line at y = graph_height
        for (int i = 0; i < WIDTH; ++i)
            image.setPixel(i, graph_height, sf::Color::White);

        double p1, p2 = 0;
        int max_height = HEIGHT - graph_height;

        // draws the contents of the lagrangians vector
        for (int index = 1; index < lagrangians.size(); index++){

            // makes sure it is a valid position
            p2 = std::clamp((int)(HEIGHT - lagrangians[index] * max_height / highest_lagrangian), 0, HEIGHT-1);
            p1 = std::clamp((int)(HEIGHT - lagrangians[index-1] * max_height / highest_lagrangian), 0, HEIGHT-1);

            //draws a pixel at x = index, y = p2
            image.setPixel(index, p2, sf::Color::White);

            int // start and stop y for the vertical line
                start = (int)std::min(p1, p2),
                stop = (int)std::max(p1, p2)           
            ;

            // connects the previous point with the current using a verical line
            // from y = start to y = stop at x = index
            for (int i = start; i < stop; i++){
                image.setPixel(index, i, sf::Color::White);
            }

        }
    }

    /*
    * draws every body in the Simulation::bodies
    * with the respective color on the "image"
    * object
    */
    void draw_minimal(){
        for (const auto& p : Simulation::bodies){
            // if it's in the canva range it draws it 
            if (p.position.x >= 0 && p.position.x < Simulation::WIDTH && p.position.y >= 0 && p.position.y < graph_height)
                image.setPixel(static_cast<unsigned int>(p.position.x), static_cast<unsigned int>(p.position.y), Body::get_color(p.temp));
        }
    }

    void draw_detailed(){
        for (const auto& p : Simulation::bodies){
            // if it's in the canva range it draws it 
            if (p.position.x >= 0 && p.position.x < WIDTH && p.position.y >= 0 && p.position.y < graph_height)
                draw_circle(std::log(p.mass)/10, p.position, p.color);
        }
        clear_graph();
        draw_graph();
    }



    render_func render_modes[] = {draw_detailed, draw_minimal, draw_minimal};

}
}