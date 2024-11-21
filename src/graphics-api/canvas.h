/*
* this file contains all the functions
* relative to image inits, drawing and
* graphic related functions 
*/

sf::Image blank;
sf::Image image;


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
void draw_circle(int radius, int centerX, int centerY){
       for (int y = centerY - radius; y <= centerY + radius; ++y) {
           for (int x = centerX - radius; x <= centerX + radius; ++x) {
               int dx = x - centerX;
               int dy = y - centerY;
               if (dx * dx + dy * dy <= radius * radius) {
                   if (x >= 0 && x < (int)image.getSize().x && y >= 0 && y < (int)image.getSize().y) {
                       image.setPixel(x, y, sf::Color::White);
                   }
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

/*
* draws every body in the Simulation::bodies
* with the respective color on the "image"
* object
*/
void draw_bodies(){
        for (const auto& p : Simulation::bodies){
            // if it's in the canva range it draws it 
            if (p.position.x >= 0 && p.position.x < Simulation::WIDTH && p.position.y >= 0 && p.position.y < Simulation::HEIGHT)
                image.setPixel(static_cast<unsigned int>(p.position.x), static_cast<unsigned int>(p.position.y), p.color);
        }
}