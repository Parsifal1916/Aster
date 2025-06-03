#include "Aster/impl/barnes_hut.tpp"
#include "Aster/impl/BHT_impl.tpp"

namespace Aster{
namespace Barnes{

double variation(){
    return ((std::rand() % 100) + 50)/150; 
}

template <>
vec3 pick_newpos(Simulation<vec3>* _s){
    return _s -> get_center() * variation();
}

template <>
vec2 pick_newpos(Simulation<vec2>* _s){
    return _s -> get_center() * variation();
}

/*                   //===---------------------------------------------------------===//
.                    // MORTON IMPLEMENTATION                                         //
.                    //===---------------------------------------------------------===*/

template <> 
void compare_bounding_vectors(const vec2& to_compare, vec2& res){
    if (std::abs(to_compare.x) > std::abs(res.x)) res.x = std::abs(to_compare.x);
    if (std::abs(to_compare.y) > std::abs(res.y)) res.y = std::abs(to_compare.y);
}

template <> 
void compare_bounding_vectors(const vec3& to_compare, vec3& res){
    if (std::abs(to_compare.x) > res.x) res.x = std::abs(to_compare.x);
    if (std::abs(to_compare.y) > res.y) res.y = std::abs(to_compare.y);
    if (std::abs(to_compare.z) > res.z) res.z = std::abs(to_compare.z);
}


template <> 
uint32_t get_morton<vec3>(Barnes_Hut<vec3>* _s, vec3 point) {
    point.x = (point.x) / (_s -> bounding_box.x);
    point.y = (point.y) / (_s -> bounding_box.y);
    point.z = (point.z) / (_s -> bounding_box.z);
    uint32_t ix = std::min((uint32_t)(point.x * (1 << PRECISION_BITS)), (1u << PRECISION_BITS) - 1);
    uint32_t iy = std::min((uint32_t)(point.y * (1 << PRECISION_BITS)), (1u << PRECISION_BITS) - 1);
    uint32_t iz = std::min((uint32_t)(point.z * (1 << PRECISION_BITS)), (1u << PRECISION_BITS) - 1);
    return (interleap_coord(iy, 2) << 2) | (interleap_coord(ix, 2) << 1) | interleap_coord(iz, 2); 
}

template <> 
uint32_t get_morton<vec2>(Barnes_Hut<vec2>* _s, vec2 point) {
    point.x = (point.x) / (_s -> bounding_box.x);
    point.y = (point.y) / (_s -> bounding_box.y);
    uint32_t ix = std::min((uint32_t)(point.x * (1 << PRECISION_BITS)), (1u << PRECISION_BITS) - 1);
    uint32_t iy = std::min((uint32_t)(point.y * (1 << PRECISION_BITS)), (1u << PRECISION_BITS) - 1);
    return (interleap_coord(iy) << 1) | interleap_coord(ix); 
}

/**
* @brief interleaps the bits of a given number
*/
uint32_t interleap_coord(uint16_t num, size_t shift){
    auto x = (uint32_t)num;
    if (shift == 1){
        x = (x | (x << 8)) & 0x00FF00FF;
        x = (x | (x << 4)) & 0x0F0F0F0F;
        x = (x | (x << 2)) & 0x33333333;
        x = (x | (x << 1)) & 0x55555555;
    } else {
        x &= 0x3FF;  
        x = (x | (x << 16)) & 0x030000FF;
        x = (x | (x << 8))  & 0x0300F00F;
        x = (x | (x << 4))  & 0x030C30C3;
        x = (x | (x << 2))  & 0x09249249;
    }
    return x;
}

}
}