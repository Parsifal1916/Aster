#include <iomanip>
#include <iostream>

#include "Aster/graphics/3d_graphics.h"
#include "Aster/graphics/inferno_scale.h"
#include "Aster/graphics/clienting.h"
//#include "graphics/2d_graphics.h"

#include "Aster/simulations/BHT_sim.h"
#include "Aster/simulations/barnes-hut3d.h"
#include "Aster/simulations/barnes-hut.h"

namespace Aster{
namespace Renderer{

//===---------------------------------------------------------===//
// 3d rendering impl                                             //
//===---------------------------------------------------------===//

vec3 rotate_point(Simulation3d* _s, vec3 v){
    v = v - rot_center;
    float x1, z1, y1, z2;

    x1 = v.x * cos_x_theta + v.z * sin_x_theta;
    z1 = - v.x * sin_x_theta + v.z *cos_x_theta;

    y1 = v.y * cos_y_theta + z1 *sin_y_theta;
    z2 = -v.y *sin_y_theta + z1 * cos_y_theta;

    float scale = distance/ (distance + z2);

    v.x = x1*scale; 
    v.y = y1*scale; 
    v.z = z2;

    v.x *= _s -> data.depth / distance;
    v.y *= _s -> data.depth / distance;

    return v + rot_center;
}

bool is_in3d_bounds(Simulation3d* _s, vec3 v){
    return v.x >= 0 && v.x < _s -> data.WIDTH && 
           v.y >= 0 && v.y < _s -> data.HEIGHT;
}

void clear_graph3d(Simulation3d* _s){
    for (int i = 0; i < _s -> data.WIDTH; i++){
    for (int j = _s -> data.graph_height; j < _s -> data.HEIGHT; j++){
        clients3d[_s].screen.setPixel(i, j, bg_color);
    }
    }
}

void draw_graph3d(Simulation3d* _s){
    std::ostringstream scien;
    scien << std::scientific << std::setprecision(3) << _s -> lagrangians.back();

    clients3d[_s].lagr_text.setString("Lagrangian: "+ scien.str() + "J");
    clients3d[_s].lagr_text.setPosition(10, _s -> data.graph_height + 10);
    clients3d[_s].window -> draw(clients3d[_s].lagr_text);

    for (int i = 0; i < _s -> data.WIDTH; ++i)
        clients3d[_s].screen.setPixel(i, _s -> data.graph_height, sf::Color::White);

    double p1, p2 = 0;
    int max_height = _s -> data.HEIGHT - _s -> data.graph_height;

    for (int index = 1; index < _s -> lagrangians.size(); index++){
        p2 = std::clamp(_s -> data.HEIGHT - _s -> lagrangians[index] * max_height / _s -> highest_lagrangian, 0.0, _s -> data.HEIGHT-1);
        p1 = std::clamp(_s -> data.HEIGHT - _s -> lagrangians[index-1] *  max_height / _s -> highest_lagrangian, 0.0, _s -> data.HEIGHT-1);
        
        //assert(p2 < HEIGHT);
        clients3d[_s].screen.setPixel(index, p2, sf::Color::White);
        
        int 
            start = (int)std::min(p1, p2),
            stop = (int)std::max(p1, p2)           
        ;
        
        for (int i = start; i < stop; i++){
            clients3d[_s].screen.setPixel(index, i, sf::Color::White);
        }
    }
}

void draw_termal3d(Simulation3d* _s){
    vec3 temp = {0,0,0};
    for (const auto& p : _s -> bodies){
        // if it's in the canva range it draws it 
        temp = rotate_point(_s, p.position);
        if (is_in3d_bounds(_s, temp))
            clients3d[_s].screen.setPixel(static_cast<unsigned int>(temp.x), static_cast<unsigned int>(temp.y), sf::Color::White);
    }
}

void draw_sphere(int radius, vec3 center, Simulation3d* _s, sf::Color color){
    for (int y = center.y - radius; y <= center.y + radius; ++y) {
        for (int x = center.x - radius; x <= center.x + radius; ++x) {
            int dx = x - center.x;
            int dy = y - center.y;
            if (dx * dx + dy * dy <= radius * radius) {
                if (x >= 0 && x < _s -> data.WIDTH && y >= 0 && y < _s -> data.HEIGHT)
                    clients3d[_s].screen.setPixel(x, y, color);
            }
        }
    }
}

void clear_graph(Simulation3d* _s){
    for (int i = 0; i < _s -> data.WIDTH; i++){
    for (int j = _s -> data.graph_height; j < _s -> data.HEIGHT; j++){
        clients3d[_s].screen.setPixel(i, j, bg_color);
    }
    }
}

void handle_displacement(){
    if (inputs["up"])
        y_theta += DIS_SCALE;
    if (inputs["down"])    
        y_theta -= DIS_SCALE;
    if (inputs["left"])
        x_theta -= DIS_SCALE;
    if (inputs["right"])
        x_theta += DIS_SCALE;

    cos_x_theta = cos(x_theta);
    cos_y_theta = cos(y_theta);

    sin_x_theta = sin(x_theta);
    sin_y_theta = sin(y_theta); 
}


void reset(Simulation3d* _s){
    clients3d[_s].screen = clients3d[_s].blank;
}

/*
* draws every body in the Simulation::bodies
* with white on the "image"
* object
*/
void draw_minimal3d(Simulation3d* _s){
    vec3 temp = {0,0,0};
    for (const auto& p : _s -> bodies){
        // if it's in the canva range it draws it 
        temp = rotate_point(_s, p.position);
        if (is_in3d_bounds(_s, temp) )
            clients3d[_s].screen.setPixel(
                static_cast<unsigned int>(temp.x), 
                static_cast<unsigned int>(temp.y), 
            inferno_cm[int(std::max(0.0, std::min(1.0, p.temp/10))*255)]
        );
    }
}

void draw_detailed3d(Simulation3d* _s){
    vec3 temp = {0,0,0};
    for (const auto& p : _s -> bodies){
        // if it's in the canva range it draws it 
        temp = rotate_point(_s, p.position);
        if (temp.x >= 0 && temp.x < _s -> data.WIDTH && temp.y >= 0 && temp.y < _s -> data.graph_height)
            draw_sphere(std::log(p.mass)/10, temp, _s, sf::Color::White);
    }
    if (!_s -> data.show_graph) return;
    clear_graph(_s);
    draw_graph3d(_s);
}

void serve_new3d(Simulation3d* _s){
    Client client;

    if (!clients3d.size()) Renderer::setup();
    sf::Image f, b; 
    client.screen.create(_s -> data.WIDTH, _s -> data.HEIGHT,  bg_color);
    client.blank.create(_s -> data.WIDTH,_s -> data.HEIGHT , bg_color);
    client.window  = new sf::RenderWindow();
    client.window -> create(sf::VideoMode(_s -> data.WIDTH, _s -> data.HEIGHT), "Aster's Simulation");
    client.is_3d =  true;

    clients3d[_s] = client;
    if (!_s -> data.show_graph) return;

    sf::Text& lagr_text = clients3d[_s].lagr_text;

    lagr_text.setFont(font);
    lagr_text.setString("Langragiano: ");
    lagr_text.setCharacterSize(16);
    lagr_text.setFillColor(sf::Color::White);
}

inline static bool is_open(Simulation3d* _s){
    return clients3d[_s].window -> isOpen();
}

using render_func3d = void(*)(Simulation3d*);
render_func3d render3d;
render_func3d render_modes3d[] = {draw_minimal3d, draw_minimal3d, draw_minimal3d, draw_minimal3d, draw_termal3d};

void handle_keyboard_input(sf::Event e){
    if (e.type == sf::Event::KeyPressed) {
        switch (e.key.code) {
            case sf::Keyboard::Up:
                inputs["up"] = true;
                break;
            case sf::Keyboard::Down:
                inputs["down"] = true;
                break;
            case sf::Keyboard::Left:
                inputs["left"] = true;
                break;
            case sf::Keyboard::Right:
                inputs["right"] = true;
                break;
            case sf::Keyboard::Space:
                inputs["space"] = true;
                paused = !paused;
                break;
            default:
                break;
        }
    }
    else if (e.type == sf::Event::KeyReleased) {
        switch (e.key.code) {
            case sf::Keyboard::Up:
                inputs["up"] = false;
                break;
            case sf::Keyboard::Down:
                inputs["down"] = false;
                break;
            case sf::Keyboard::Left:
                inputs["left"] = false;
                break;
            case sf::Keyboard::Right:
                inputs["right"] = false;
                break;
            case sf::Keyboard::Space:
                inputs["space"] = false;
                break;
            default:
                break;
        }
    }
    else if (e.type == sf::Event::MouseWheelScrolled){
        if (e.mouseWheelScroll.wheel != sf::Mouse::VerticalWheel) return;
        distance += e.mouseWheelScroll.delta*100;
    }
}

void do_window3d(Simulation3d* _s){
    sf::Texture texture;
    sf::Sprite sprite;

    Renderer::serve_new3d(_s);
    sf::RenderWindow& window = *Renderer::clients3d[_s].window;
    distance = _s -> data.depth*2;
    
    std::cout << "\n";
    
    auto tot_time = std::chrono::high_resolution_clock::now()-std::chrono::high_resolution_clock::now();
    auto sec_time = std::chrono::high_resolution_clock::now()-std::chrono::high_resolution_clock::now();
    
    while (is_open(_s)){
        auto start = std::chrono::high_resolution_clock::now();
        sf::Event event;
    
        while (window.pollEvent(event)){
            if (event.type == sf::Event::Closed) window.close();
            handle_keyboard_input(event);
        }    
    
        handle_displacement();            
        render_modes3d[_s -> data.type](_s);
        texture.loadFromImage(Renderer::clients3d[_s].screen);
        sprite.setTexture(texture, true);   
        window.clear();
        window.draw(sprite);
        window.display();
    
        if (!_s -> data.leave_traces)
            Renderer::reset(_s);
        
        if (_s -> data.save)
            texture.copyToImage().saveToFile("saves/file" + std::to_string((int)(_s -> time_passed)) + ".png");
        
        auto start_second = std::chrono::high_resolution_clock::now();
    
        if (!paused) {
        _s -> step();
        _s -> time_passed++;              
        }
    
        rot_center = vec3({_s -> data.WIDTH / 2, _s -> data.HEIGHT /2,_s -> data.depth/2});

        auto end = std::chrono::high_resolution_clock::now();
        tot_time += end-start;
        sec_time += end-start_second;
    
        std::cout << "Totale: " 
            << std::chrono::duration_cast<std::chrono::milliseconds>(tot_time).count() / _s -> time_passed << "ns\nSimulazione: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(sec_time).count() / _s -> time_passed << "           "
        ; 
    }
}

void save_stepbystep(Simulation3d* _s){
    sf::Texture texture;
    sf::Sprite sprite;

    serve_new3d(_s);
        
    //sf::RenderTexture render_txture;
    //render_txture.create(_s -> data.WIDTH, _s -> data.HEIGHT);
    sf::RenderWindow& window = *clients3d[_s].window;
    while (is_open(_s)){
        sf::Event event;
        while (window.pollEvent(event)){
            if (event.type == sf::Event::Closed) window.close();
        }    

        render_modes3d[_s -> data.type](_s);
        texture.loadFromImage(clients3d[_s].screen);

        sprite.setTexture(texture, true);
        texture.copyToImage().saveToFile("file" + std::to_string((int)(_s -> time_passed)) + ".png");
        
        window.clear();
        window.draw(sprite);
        window.display();

        if (!_s -> data.leave_traces)
            Renderer::reset(_s);
        
        _s -> step();
        _s -> time_passed++;
    }
}

//===---------------------------------------------------------===//
// 2d rendering impl                                             //
//===---------------------------------------------------------===//

std::unordered_map<Simulation*, Client> clients;

void draw_circle(int radius, vec2 center, Simulation* _s, sf::Color color = sf::Color::White){
    for (int y = center.y - radius; y <= center.y + radius; ++y) {
        for (int x = center.x - radius; x <= center.x + radius; ++x) {
            int dx = x - center.x;
            int dy = y - center.y;
            if (dx * dx + dy * dy <= radius * radius) {
                if (x >= 0 && x < _s -> data.WIDTH && y >= 0 && y < _s -> data.HEIGHT)
                    clients[_s].screen.setPixel(x, y, color);
            }
        }
    }
}

void clear_graph(Simulation* _s){
    for (int i = 0; i < _s -> data.WIDTH; i++){
    for (int j = _s -> data.graph_height; j < _s -> data.HEIGHT; j++){
        clients[_s].screen.setPixel(i, j, bg_color);
    }
    }
}

void draw_graph(Simulation* _s){
    std::ostringstream scien;
    scien << std::scientific << std::setprecision(3) << _s -> lagrangians.back();
    clients[_s].lagr_text.setString("Lagrangian: "+ scien.str() + "J");
    clients[_s].lagr_text.setPosition(10, _s -> data.graph_height + 10);
    clients[_s].window -> draw(clients[_s].lagr_text);

    for (int i = 0; i < _s -> data.WIDTH; ++i)
        clients[_s].screen.setPixel(i, _s -> data.graph_height, sf::Color::White);

    double p1, p2 = 0;
    int max_height = _s -> data.HEIGHT - _s -> data.graph_height;

    for (int index = 1; index < _s -> lagrangians.size(); index++){
        p2 = std::clamp(_s -> data.HEIGHT - _s -> lagrangians[index] * max_height / _s -> highest_lagrangian, 0.0, _s -> data.HEIGHT-1);
        p1 = std::clamp(_s -> data.HEIGHT - _s -> lagrangians[index-1] *  max_height / _s -> highest_lagrangian, 0.0, _s -> data.HEIGHT-1);

        //assert(p2 < HEIGHT);
        clients[_s].screen.setPixel(index, p2, sf::Color::White);

        int 
            start = (int)std::min(p1, p2),
            stop = (int)std::max(p1, p2)           
        ;

        for (int i = start; i < stop; i++){
            clients[_s].screen.setPixel(index, i, sf::Color::White);
        }

    }
}

void reset(Simulation* _s){
    clients[_s].screen = clients[_s].blank;
}

int clamp_rgb(double x){
    return (int)std::max(std::min(int(x), 255), 0);
}

void draw_termal(Simulation* _s){
    for (const auto& p : _s -> bodies){
        // if it's in the canva range it draws it  
        if (p.position.x >= 0 && p.position.x < _s -> data.WIDTH && p.position.y >= 0 && p.position.y < _s -> data.HEIGHT)
            clients[_s].screen.setPixel(static_cast<unsigned int>(p.position.x), static_cast<unsigned int>(p.position.y), sf::Color::White);
    }
}

void draw_minimal(Simulation* _s){
    for (const auto& p : _s -> bodies){
        // if it's in the canva range it draws it 
        if (p.position.x >= 0 && p.position.x < _s -> data.WIDTH && p.position.y >= 0 && p.position.y < _s -> data.HEIGHT)
            clients[_s].screen.setPixel(static_cast<unsigned int>(p.position.x), static_cast<unsigned int>(p.position.y), sf::Color::White);
    }
}

void draw_detailed(Simulation* _s){
    for (const auto& p : _s -> bodies){
        
        // if it's in the canva range it draws it 
        if (p.position.x >= 0 && p.position.x < _s -> data.WIDTH && p.position.y >= 0 && p.position.y < _s -> data.graph_height)
            draw_circle(std::log(p.mass)/10, p.position, _s, sf::Color::White);
    }
    if (!_s -> data.show_graph) return;
    clear_graph(_s);
    draw_graph(_s);
}

using render_func = void(*)(Simulation*);
render_func render;
render_func render_modes[] = {draw_minimal, draw_minimal, draw_minimal, draw_minimal, draw_termal};

void serve_new(Simulation* _s){
    Client client;
    
    if (!clients.size()) Renderer::setup();
    sf::Image f, b; 

    client.screen.create(_s -> data.WIDTH, _s -> data.HEIGHT,  bg_color);
    client.blank.create(_s -> data.WIDTH,_s -> data.HEIGHT , bg_color);
    client.window  = new sf::RenderWindow();
    client.window -> create(sf::VideoMode(_s -> data.WIDTH, _s -> data.HEIGHT), "Aster's Simulation");
    clients[_s] = client;
    
    if (!_s -> data.show_graph) return;
    sf::Text& lagr_text = clients[_s].lagr_text;
    
    lagr_text.setFont(font);
    lagr_text.setString("Langragiano: ");
    lagr_text.setCharacterSize(16);
    lagr_text.setFillColor(sf::Color::White);
}

inline static bool is_open(Simulation* _s){
    return clients[_s].window -> isOpen();
}

void do_window(Simulation* _s){
    sf::Texture texture;
    sf::Sprite sprite;

    Renderer::serve_new(_s);
        
    //sf::RenderTexture render_txture;
    //render_txture.create(_s -> data.WIDTH, _s -> data.HEIGHT);

    sf::RenderWindow& window = *Renderer::clients[_s].window;

    while (Renderer::is_open(_s)){
        sf::Event event;
        while (window.pollEvent(event)){
            if (event.type == sf::Event::Closed) window.close();
        }    

        Renderer::render_modes[_s -> data.type](_s);

        texture.loadFromImage(Renderer::clients[_s].screen);
        sprite.setTexture(texture, true);

        if (_s -> data.save)
            texture.copyToImage().saveToFile("file" + std::to_string((int)(_s -> time_passed)) + ".png");

        window.clear();
        window.draw(sprite);
        window.display();

        if (!_s -> data.leave_traces)
            Renderer::reset(_s);
        
        _s -> step();
        _s -> time_passed++;
    }
}

void save_stepbystep(Simulation* _s){
    sf::Texture texture;
    sf::Sprite sprite;

    Renderer::serve_new(_s);
        
    //sf::RenderTexture render_txture;
    //render_txture.create(_s -> data.WIDTH, _s -> data.HEIGHT);

    sf::RenderWindow& window = *Renderer::clients[_s].window;

    while (Renderer::is_open(_s)){
        sf::Event event;
        while (window.pollEvent(event)){
            if (event.type == sf::Event::Closed) window.close();
        }    

        Renderer::render_modes[_s -> data.type](_s);

        texture.loadFromImage(Renderer::clients[_s].screen);
        sprite.setTexture(texture, true);

        texture.copyToImage().saveToFile("file" + std::to_string((int)(_s -> time_passed)) + ".png");

        window.clear();
        window.draw(sprite);
        window.display();

        if (!_s -> data.leave_traces)
            Renderer::reset(_s);
        
        _s -> step();
        _s -> time_passed++;
    }
}

}
}