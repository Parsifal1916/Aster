#include <cmath>

#include "Aster/physics/tool-chain.h"
#include "Aster/physics/vectors.h"

#include "Aster/simulations/sim_obj.h"
#include "Aster/simulations/3d_sim_obj.h"

namespace Aster{
//===---------------------------------------------------------===//
// Update methods                                                //
//===---------------------------------------------------------===//

void update_euler(Body* b, Simulation* _s){
    b -> velocity += b -> acceleration* _s -> data.dt;
    b -> position += b -> velocity * _s -> data.dt;
    b -> acceleration = {0, 0};
}

void update_leapfrog(Body* b, Simulation* _s){
   b -> position += b -> velocity * _s -> data.dt + b -> acceleration * _s -> data.dt * _s -> data.dt * .5;
   b -> velocity += (b -> acceleration + b -> prev_acc)* _s -> data.dt * .5;
   b -> acceleration =  {0,0};
}

void update_symplectic4(Body* body, Simulation* _s){
    body -> position += body -> velocity * c1 * _s -> data.dt;
    
    _s -> update_pair(body);
    body -> velocity += body -> acceleration * d1 *_s -> data.dt;
    body -> position += body -> velocity * c2 *_s -> data.dt;
    
    _s -> update_pair(body);
    body -> velocity += body -> acceleration * d2 *_s -> data.dt;
    body -> position += body -> velocity * c2 *_s -> data.dt;
    
    _s -> update_pair(body);
    body -> velocity += body -> acceleration * d1 *_s -> data.dt;
    body -> position += body -> velocity * c1 * _s -> data.dt;
}

//===---------------------------------------------------------===//
// Force calculation methods                                     //
//===---------------------------------------------------------===//

vec2 newtonian(double m1, double m2, vec2 v1, vec2 v2, vec2 p1, vec2 p2, Simulation* _s){
    vec2 d = p2 - p1;
    return d.normalize() *_s -> data.G* m1*m1/ (d.sqr_magn()+ 1);
}

vec2 pn2(double m1, double m2, vec2 v1, vec2 v2, vec2 p1, vec2 p2, Simulation* _s){ 
    vec2 n = p2 - p1;
    double r = n.magnitude() * _s -> data.simulation_scale+1;

    double term1 = _s -> data.G*m2/(_s -> data.c_squared*r*r);
    double m_ratio = m1/m2;
    vec2 v = v1 - v2;

    double correction = std::exp(-r*r/_s -> data.takeover)*_s -> data.e_squared;
    vec2 a_newton = n * _s -> data.G *m1*m2/(r*r);
    vec2 pn1 = ( n*(4 + 2*m_ratio)*term1*_s -> data.c_squared - v*n*(1+ 3*m_ratio) + v*n*v*4 - v*n*n*3/2) * term1 ;
    vec2 pn25 = (v*n)*v*-8/5 * _s -> data.G* _s -> data.G*m1*m2*(m1+m2)/(_s -> data.c_squared * _s -> data.c_squared* _s -> data.c*r*r*r);

    return a_newton + pn1;// + pn25 + n*correction;
}

std::vector<vec2> get_new_pos(vec2* position, vec2* velocity, vec2* acceleration, double step){
    std::vector<vec2> retval;
    retval.push_back(vec2({velocity -> x + acceleration -> x *step, velocity -> y + acceleration -> y *step}));
    retval.push_back(vec2({position -> x + retval[0].x * step     , position -> y + retval[0].y * step     }));
    return retval;
}
vec2 rk4(double m1, double m2, vec2 v1, vec2 v2, vec2 p1, vec2 p2, Simulation* _s){
    vec2 k1, k2, k3, k4;
    std::vector<vec2> dummy;
    k1 = _s -> get_force(m1, m2, v1, v2, p1, p2, _s) * _s -> data.dt;
    dummy = get_new_pos(&p1, &v1, &k1, _s -> data.dt/2);
    k2 = _s -> get_force(m1, m2, v1 + dummy[0]/2.f, v2, p1 + dummy[1]/2.f, p2, _s) * _s -> data.dt;
    dummy = get_new_pos(&p1, &v1, &k2, _s -> data.dt/2);
    k3 = _s -> get_force(m1, m2, v1 + dummy[0]/2.f, v2, p1 + dummy[1]/2.f, p2, _s) * _s -> data.dt;
    dummy = get_new_pos(&p1, &v1, &k3, _s -> data.dt);
    k4 = _s -> get_force(m1, m2, v1 + dummy[0]/2.f, v2, p1 + dummy[1]/2.f, p2, _s) * _s -> data.dt;

    return (k1 + k2*2.f + k3*2.f + k4)/6.f;
}

//===---------------------------------------------------------===//
// Update methods (3d)                                           //
//===---------------------------------------------------------===//


void update_euler_3d(Body3d* b, Simulation3d* _s){
    b -> velocity += b -> acceleration* _s -> data.dt;
    b -> position += b -> velocity * _s -> data.dt;
    b -> acceleration = {0, 0, 0};
}

void update_leapfrog_3d(Body3d* b, Simulation3d* _s){
   b -> position += b -> velocity * _s -> data.dt + b -> acceleration * _s -> data.dt * _s -> data.dt * .5;
   b -> velocity += (b -> acceleration + b -> prev_acc)* _s -> data.dt * .5;
   b -> acceleration =  {0, 0, 0};
}

void update_symplectic4_3d(Body3d* body,Simulation3d* _s){
    body -> position += body -> velocity * c1 * _s -> data.dt;
    
    _s -> update_pair(body);
    body -> velocity += body -> acceleration * d1 *_s -> data.dt;
    body -> position += body -> velocity * c2 *_s -> data.dt;
    
    _s -> update_pair(body);
    body -> velocity += body -> acceleration * d2 *_s -> data.dt;
    body -> position += body -> velocity * c2 *_s -> data.dt;
    
    _s -> update_pair(body);
    body -> velocity += body -> acceleration * d1 *_s -> data.dt;
    body -> position += body -> velocity * c1 * _s -> data.dt;
}

//===---------------------------------------------------------===//
// Force calculation methods (3d)                                //
//===---------------------------------------------------------===//

vec3 newtonian_3d(double m1, double m2, vec3 v1, vec3 v2, vec3 p1, vec3 p2, Simulation3d* _s){
    vec3 d = p2 - p1;
    return d.normalize() *_s -> data.G* m1*m1/ (d.sqr_magn() + 1);
}

/*
* evaluates the post-newtonian approx.
* @param m1: mass of the first object
* @param m2: mass of the second object
* @param v1: velocity of the first object
* @param v2: velocity of the second object 
* @param p1: position of the first object
* @param p2: position if the second object
*/
vec3 pn2_3d(double m1, double m2, vec3 v1, vec3 v2, vec3 p1, vec3 p2, Simulation3d* _s){ 
    vec3 n = p2 - p1;
    double r = n.magnitude() * _s -> data.simulation_scale+1;

    double term1 = _s -> data.G*m2/(_s -> data.c_squared*r*r);
    double m_ratio = m1/m2;
    vec3 v = v1 - v2;

    double correction = std::exp(-r*r/_s -> data.takeover)*_s -> data.e_squared;
    vec3 a_newton = n * _s -> data.G *m1*m2/(r*r);
    vec3 pn1 = ( n*(4 + 2*m_ratio)*term1*_s -> data.c_squared - v*n*(1+ 3*m_ratio) + v*n*v*4 - v*n*n*3/2) * term1 ;
    vec3 pn25 = (v*n)*v*-8/5 * _s -> data.G* _s -> data.G*m1*m2*(m1+m2)/(_s -> data.c_squared * _s -> data.c_squared* _s -> data.c*r*r*r);

    return a_newton + pn1;// + pn25 + n*correction;
}

std::vector<vec3> get_new_pos_3d(vec3* position, vec3* velocity, vec3* acceleration, double step){
    std::vector<vec3> retval;
    retval.push_back(velocity -> update_by(acceleration, step));
    retval.push_back(position -> update_by(&retval[0], step));
    return retval;
}
vec3 rk4_3d(double m1, double m2, vec3 v1, vec3 v2, vec3 p1, vec3 p2, Simulation3d* _s){
    vec3 k1, k2, k3, k4;
    std::vector<vec3> dummy;
    k1 = _s -> get_force(m1, m2, v1, v2, p1, p2, _s) * _s -> data.dt;
    dummy = get_new_pos_3d(&p1, &v1, &k1, _s -> data.dt/2);
    k2 = _s -> get_force(m1, m2, v1 + dummy[0]/2.f, v2, p1 + dummy[1]/2.f, p2, _s) * _s -> data.dt;
    dummy = get_new_pos_3d(&p1, &v1, &k2, _s -> data.dt/2);
    k3 = _s -> get_force(m1, m2, v1 + dummy[0]/2.f, v2, p1 + dummy[1]/2.f, p2, _s) * _s -> data.dt;
    dummy = get_new_pos_3d(&p1, &v1, &k3, _s -> data.dt);
    k4 = _s -> get_force(m1, m2, v1 + dummy[0]/2.f, v2, p1 + dummy[1]/2.f, p2, _s) * _s -> data.dt;

    return (k1 + k2*2.f + k3*2.f + k4)/6.f;
}

//===---------------------------------------------------------===//
// Cosmology stuff                                               //
//===---------------------------------------------------------===//

/*
* gets the delta a from the previous 
* value of the scale factor 
*/
float get_da(float s, Simulation* _s){
    double term1 = _s -> data.cosmological_constant * s *s / 3 + 10e-11;

    double term2 = _s -> data.field_eq_scalar * _s -> data.c_squared * s *s /3 * _s -> obj / (_s -> data.HEIGHT * _s -> data.WIDTH * s * s * _s -> data.simulation_scale * _s -> data.simulation_scale);
    return _s -> data.c * std::sqrt(term2  + term1 );
}

/*
* updates the scale factor with the
* rk4 method
*/
void update_scale(Simulation* _s){
    double k1 = get_da(_s -> data.current_a, _s) * _s -> data.dt;
    double k2 = get_da(_s -> data.current_a + k1 * 0.5, _s) * _s -> data.dt;
    double k3 = get_da(_s -> data.current_a + k2 * 0.5, _s) * _s -> data.dt;
    double k4 = get_da(_s -> data.current_a + k3, _s) * _s -> data.dt;

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


vec2 vec2::operator*(const vec2& other) const{
    return vec2(x * other.x, y * other.y);
}

vec2& vec2::operator/=(const double& other){
    this -> x /= other;
    this -> y /= other;
    return *this;
}

float vec2::magnitude(){
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

double vec2::sqr_magn() const {
    return x*x + y*y;
}

// vec3

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

vec3 vec3::operator*(const vec3& other) const{
    return vec3(x * other.x, y * other.y, z*other.z);
}

vec3& vec3::operator/=(const double& other){
    this -> x /= other;
    this -> y /= other;
    this -> z /= other;
    return *this;
}

float vec3::magnitude(){
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

double vec3::sqr_magn() const{
    return x*x + y*y + z*z;
}

vec3 vec3::update_by(vec3* v, double delta){
    return {
        x + v -> x *delta,
        y + v -> y *delta,
        z + v -> z *delta
    };
}


}