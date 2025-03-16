#pragma once 

#include <GLFW/glfw3.h>
#include <iostream>
#include <vector>

namespace Aster{
namespace Text{

    /**
    * @brief alerts nano vg of a resive event
    */
    void nvg_resize(int w, int h);

    /**
    * @brief ends nvg frame
    */
    void end_vg_frame();

    /**
    * @brief loads the font texture, called at the start of the renderer
    */
    void load_font();

    /**
    * @brief renders some text on screen
    * @param x: x coordinate of the text
    * @param y: y coordinate of the text
    * @param size: font size
    * @param text: text to render
    */
    void write2screen(float x, float y, float size, const char* text);

}
}