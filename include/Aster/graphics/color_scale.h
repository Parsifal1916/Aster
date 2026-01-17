#pragma once
#include <vector>

namespace Aster{
namespace Renderer{

extern float color_scale[256][3];
extern float jet_color_scale[256][3];

double get_coloring_index(double temp);


}
}