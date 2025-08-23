#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "Aster/graphics/text_rendering.h"
#include "Aster/building-api/logging.h"

#define GLFW_INCLUDE_GLEXT
#define NANOVG_GL3_IMPLEMENTATION
#include "Aster/thirdparty/nanovg.h"
#include "Aster/thirdparty/nanovg_gl.h"

#include <iostream>
#include <map>
#include <cstdlib> 
#include <vector>

namespace Aster{
namespace Text{

NVGcontext* vg;

/**
* @brief loads the font texture, called at the start of the renderer
*/
void load_font() {
    // better text quality and anti-aliasing
    glfwWindowHint(GLFW_SAMPLES, 4);
    glEnable(GL_MULTISAMPLE);

    // laods vg and checks for errors
    vg = nvgCreateGL3(NVG_ANTIALIAS | NVG_STENCIL_STROKES); 
    if (critical_if(vg == NULL, "could not init NanoVg"))
        exit(-1);


    // loads the font and checks for errors
    int font = nvgCreateFont(vg, "sans", (std::string(LIBRARY_PATH) + "/roboto.ttf").c_str());
    if (critical_if(font == -1, "could not load font, you might need to reinstall the library"))
        exit(-1);
}

/**
* @brief alerts nano vg of a resive event
*/
void nvg_resize(int w, int h){
    nvgBeginFrame(vg, w, h, 2.0f);
}

/**
* @brief ends nvg frame
*/
void end_vg_frame(){
    nvgEndFrame(vg);
}

/**
* @brief renders some text on screen
* @param x: x coordinate of the text (normalized)
* @param y: y coordinate of the text (normalized)
* @param size: font size
* @param text: text to render
*/
void write2screen(float x, float y, float scale, const char* text) {
    // sets scale and type of font
    nvgFontSize(vg, scale);
    nvgFontFace(vg, "sans");

    // sets alignment and color
    nvgFillColor(vg, nvgRGBA(255, 255, 255, 255));
    nvgTextAlign(vg, NVG_ALIGN_LEFT | NVG_ALIGN_MIDDLE);

    // shows the text at the given position
    nvgText(vg, x, y, text, NULL);
}

}
}