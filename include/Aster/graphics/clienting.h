#pragma once

namespace Aster{
namespace Renderer{


sf::Font font;

struct Client{
    sf::Image blank, screen;
    sf::RenderWindow* window;
    sf::Text lagr_text;
    bool display_graph = false;
    bool is_3d = false;

     Client() = default;
};

std::vector<sf::Color> rng_colors = {
    sf::Color::White,
    sf::Color::Red,
    sf::Color::Blue,
    sf::Color::Green,
};

sf::Color bg_color = sf::Color::Black;

void setup(){
    font.loadFromFile("/usr/share/fonts/TTF/Hack-Regular.ttf");
}


}
}