#include <cmath>

#include "Aster/physics/tool-chain.h"
#include "Aster/physics/vectors.h"

#include "Aster/simulations/sim_obj.h"

#include "Aster/impl/tool_chain_impl.tpp"

namespace Aster{
const REAL PI = 3.141592653589793238462643383279;

//===---------------------------------------------------------===//
// Cosmology stuff                                               //
//===---------------------------------------------------------===//

/*
* gets the delta a from the previous 
* value of the scale factor 
*/
template <typename T>
float get_da(float s, Simulation<T>* _s){
    REAL term1 = _s -> data.cosmological_constant * s *s / 3 + 10e-11;

    REAL term2 = _s -> data.field_eq_scalar * _s -> data.c_squared * s *s /3 * _s -> obj / (_s -> data.HEIGHT * _s -> data.WIDTH * s * s * _s -> data.simulation_scale * _s -> data.simulation_scale);
    return _s -> data.c * std::sqrt(term2  + term1 );
}

/*
* updates the scale factor with the
* rk4 method
*/
template <typename T> 
void update_scale(Simulation<T>* _s){
    REAL k1 = get_da<T>(_s -> data.current_a, _s) * _s -> data.dt;
    REAL k2 = get_da<T>(_s -> data.current_a + k1 * 0.5, _s) * _s -> data.dt;
    REAL k3 = get_da<T>(_s -> data.current_a + k2 * 0.5, _s) * _s -> data.dt;
    REAL k4 = get_da<T>(_s -> data.current_a + k3, _s) * _s -> data.dt;

    _s -> data.current_a += (k1 + 2 * k2 + 2 * k3 + k4) / 6;
}

//===---------------------------------------------------------===//
// Vectors                                                       //
//===---------------------------------------------------------===//

vec2 vec2::operator*(double scalar) const {
    return vec2(x * scalar, y * scalar);
}

vec2 vec2::operator/(double scalar) const {
    return vec2(x / scalar, y / scalar);
}

vec2 vec2::operator/(float scalar) const {
    return vec2(x / scalar, y / scalar);
}

vec2 vec2::operator+(const vec2& other) const {
    return vec2(x + other.x, y + other.y);
}

vec2 vec2::operator-(const vec2& other) const {
    return vec2(x - other.x, y - other.y);
}

vec2 vec2::operator-() const {
    return vec2(-x, -y);
}

void vec2::reset() {
    x = 0;
    y = 0;
}


vec2& vec2::operator+=(const vec2& other) {
    x += other.x;
    y += other.y;
    return *this;
}

bool vec2::operator==(const vec2& other) const {
    return other.x == x && other.y == y;
}

vec2 vec2::operator*(const int& other) const{
    return vec2(x*other, y*other);
}

vec2 vec2::operator/(const int& other) const{
    return vec2(x/other, y/other);
}


REAL vec2::operator*(const vec2& other) const{
    return x * other.x + y + other.y;
}

vec2& vec2::operator/=(const REAL& other){
    this -> x /= other;
    this -> y /= other;
    return *this;
}

REAL vec2::magnitude() const{
    return std::sqrt(
        this -> x * this -> x + 
        this -> y * this -> y
    );
}

vec2 vec2::direction(vec2 v2){
    return (*this)*-1 + v2;
}

vec2 vec2::normalize(){
    float magn = this -> magnitude();
    return vec2(this -> x/magn, this -> y/magn);
}

REAL vec2::sqr_magn() const {
    return x*x + y*y;
}

bool vec2::is_fine () const {
    return !std::isnan(x) && !std::isnan(x);
}

vec2 vec2::update_by(vec2* v, REAL delta){
    return {
        x + v -> x *delta,
        y + v -> y *delta
    };
}

REAL& vec2::operator[](size_t index){
    assert(index < 2);
    return index == 0 ? x : y;
}


// vec3

void vec3::reset() {
    x = 0;
    y = 0;
    z = 0;
}

REAL& vec3::operator[](size_t index){
    assert(index < 3);
    if (index == 0) return x;
    if (index == 1) return y;
    return z;
}


vec3 vec3::operator*(double scalar) const {
    return vec3(x * scalar, y * scalar, z * scalar);
}

vec3 vec3::operator/(double scalar) const {
    return vec3(x / scalar, y / scalar, z / scalar);
}

vec3 vec3::operator/(float scalar) const {
    return vec3(x / scalar, y / scalar, z / scalar);
}

vec3 vec3::operator+(const vec3& other) const {
    return vec3(x + other.x, y + other.y, z + other.z);
}

vec3 vec3::operator-(const vec3& other) const {
    return vec3(x - other.x, y - other.y, z - other.z);
}

vec3 vec3::operator-() const {
    return vec3(-x, -y, -z);
}

vec3& vec3::operator+=(const vec3& other) {
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
}

bool vec3::operator==(const vec3& other) const {
    return other.x == x && other.y == y && other.z == z;
}

vec3 vec3::operator*(const int& other) const{
    return vec3(x*other, y*other, z*other);
}

vec3 vec3::operator/(const int& other) const{
    return vec3(x/other, y/other, z/other);
} 

REAL vec3::operator*(const vec3& other) const{
    return x * other.x + y * other.y + z*other.z;
}

vec3& vec3::operator/=(const REAL& other){
    this -> x /= other;
    this -> y /= other;
    this -> z /= other;
    return *this;
}

REAL vec3::magnitude() const{
    return std::sqrt(
        this -> x * this -> x + 
        this -> y * this -> y + 
        this -> z * this -> z
    );
}

vec3 vec3::direction(vec3 v2){
    return (*this)*-1 + v2;
}

vec3 vec3::normalize(){
       float magn = this -> magnitude();
       return vec3(this -> x/magn, this -> y/magn, this -> z/magn);
}

REAL vec3::sqr_magn() const{
    return x*x + y*y + z*z;
}

vec3 vec3::update_by(vec3* v, REAL delta){
    return {
        x + v -> x *delta,
        y + v -> y *delta,
        z + v -> z *delta
    };
}


vec2 operator*(REAL scalar, vec2 v){
    return v * scalar;
}

vec3 operator*(REAL scalar, vec3 v){
    return v * scalar;
}


bool vec3::is_fine() const {
    return !std::isnan(x) && !std::isnan(x) && !std::isnan(z);
}



}