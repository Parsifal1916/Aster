#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "Aster/graphics/text_rendering.h"
#include "Aster/building-api/logging.h"

#define NANOVG_GL3_IMPLEMENTATION
#include "Aster/thirdparty/nanovg.h"
#include "Aster/thirdparty/nanovg_gl.h"

#include <string>

namespace Aster {
namespace Text {

static NVGcontext* vg = nullptr;

int get_image(const char* path, int index){
    return nvgCreateImage(vg, path, index);
}



void load_png(const char *path, int img, int height, int width) {
    int w, h; nvgImageSize(Text::vg, img, &w, &h);
    int scale =  20;

    float aspect = float(w) / float(h);

    w = float(width) / scale * aspect;
    h = float(width) / scale;

    NVGpaint p = nvgImagePattern(Text::vg, 0, height - h, w, h, 0.0f, img, 1.0f);
    nvgBeginPath(Text::vg);
    nvgRect(Text::vg, 0, height - h, w, h);
    nvgFillPaint(Text::vg, p);
    nvgFill(Text::vg);
}


void load_font() {
    glfwWindowHint(GLFW_SAMPLES, 4);
    glEnable(GL_MULTISAMPLE);

    vg = nvgCreateGL3(NVG_ANTIALIAS | NVG_STENCIL_STROKES);
    if (critical_if(vg == nullptr, "could not init NanoVG"))
        exit(-1);

    int font = nvgCreateFont(vg, "sans", (std::string(LIBRARY_PATH) + "/roboto.ttf").c_str());
    if (critical_if(font == -1, "could not load font, you might need to reinstall the library"))
        exit(-1);
}

void begin_vg_frame(int w, int h) {
    if (!vg) return;
    nvgBeginFrame(vg, w, h, 1.0f);
}

void end_vg_frame() {
    if (!vg) return;
    nvgEndFrame(vg);
}

void write2screen(float x, float y, float scale, const char* text) {
    if (!vg || !text) return;

    nvgFontSize(vg, scale);
    nvgFontFace(vg, "sans");
    nvgFillColor(vg, nvgRGBA(255, 255, 255, 255));
    nvgTextAlign(vg, NVG_ALIGN_LEFT | NVG_ALIGN_MIDDLE);
    nvgText(vg, x, y, text, nullptr);
}

void nvg_resize(int w, int h) {
    glViewport(0, 0, w, h);
}

}
}
