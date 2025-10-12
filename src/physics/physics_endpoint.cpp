#include <cmath>
#include <cassert>

#include "Aster/physics/body.h"
#include "Aster/simulations/sim_obj.h"

namespace Aster{
extern const REAL PI = 3.141592653589793238462643383279;

vec3& BodyArray::get_position_of(size_t i){
       return positions[i];
}

vec3& BodyArray::get_velocity_of(size_t i){
       return velocities[i];
}


REAL& BodyArray::get_mass_of(size_t i){
       return masses[i];
}


REAL& BodyArray::get_temp_of(size_t i){
       return temps[i];
}


vec3& BodyArray::get_acc_of(size_t i){
       return accs[i];
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

std::ostream& operator<<(std::ostream& os, const vec3& v) {
    os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
    return os;
}

std::ostream& operator<<(std::ostream& os, const vec2& v) {
    os << "(" << v.x << ", " << v.y << ")";
    return os;
}

}