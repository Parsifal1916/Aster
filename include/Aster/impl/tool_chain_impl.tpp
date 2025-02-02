#pragma once

#include <cmath>

#include "Aster/physics/tool-chain.h"
#include "Aster/physics/vectors.h"

#include "Aster/simulations/sim_obj.h"

namespace Aster{


double c1 = 1 / (2 * (2 - std::pow(2.0, 1.0/3.0)));
double c2 = (1 - pow(2.0, 1.0/3.0)) * c1;
double d1 = 1 / (2 - std::pow(2.0, 1.0/3.0));
double d2 = -std::pow(2.0, 1.0/3.0) *d1; 
//===---------------------------------------------------------===//
// Update methods                                                //
//===---------------------------------------------------------===//


template <>
void update_euler(Body<vec2>* b, Simulation<vec2>* _s){
    b -> velocity += b -> acceleration* _s -> get_dt() ;
    b -> position += b -> velocity * _s -> get_dt() ;
}

template <>
void update_leapfrog(Body<vec2>* b, Simulation<vec2>* _s){
   b -> position += b -> velocity * _s -> get_dt()  + b -> acceleration * _s -> get_dt()  * _s -> get_dt()  * .5;
   b -> velocity += (b -> acceleration + b -> prev_acc)* _s -> get_dt()  * .5;
}

template <>
void update_symplectic4(Body<vec2>* body, Simulation<vec2>* _s){
    body -> position += body -> velocity * c1 * _s -> get_dt() ;
    
    _s -> update_pair(body);
    body -> velocity += body -> acceleration * d1 *_s -> get_dt();
    body -> position += body -> velocity * c2 *_s -> get_dt();
    
    _s -> update_pair(body);
    body -> velocity += body -> acceleration * d2 *_s -> get_dt();
    body -> position += body -> velocity * c2 *_s -> get_dt();
    
    _s -> update_pair(body);
    body -> velocity += body -> acceleration * d1 *_s -> get_dt();
    body -> position += body -> velocity * c1 * _s -> get_dt() ;
}

//===---------------------------------------------------------===//
// Force calculation methods                                     //
//===---------------------------------------------------------===//

template <>
vec2 newtonian(double m1, double m2, vec2 v1, vec2 v2, vec2 p1, vec2 p2, Simulation<vec2>* _s){
    vec2 d = p2 - p1;
    return d.normalize() *_s -> get_G()* m1*m1/ (d.sqr_magn()+ 1);
}

template <>
vec2 pn2(double m1, double m2, vec2 v1, vec2 v2, vec2 p1, vec2 p2, Simulation<vec2>* _s){ 
    vec2 n = p2 - p1;
    double r = n.magnitude() * _s -> get_scale()+1;

    double term1 = _s -> get_G()*m2/(_s -> get_c_sqr()*r*r);
    double m_ratio = m1/m2;
    vec2 v = v1 - v2;

    double correction = std::exp(-r*r/_s -> get_takeover())*_s -> get_e_sqr();
    vec2 a_newton = n * _s -> get_G() *m1*m2/(r*r);
    vec2 pn1 = ( n*(4 + 2*m_ratio)*term1*_s -> get_c_sqr() - v*n*(1+ 3*m_ratio) + v*n*v*4 - v*n*n*3/2) * term1 ;
    vec2 pn25 = (v*n)*v*-8/5 * _s -> get_G()* _s -> get_G()*m1*m2*(m1+m2)/(_s -> get_c_sqr() * _s -> get_c_sqr()* _s -> get_c()*r*r*r);

    return a_newton + pn1 + pn25 + n*correction;
}

template <>
std::vector<vec2> get_new_pos(vec2* position, vec2* velocity, vec2* acceleration, double step){
    std::vector<vec2> retval;
    retval.push_back(velocity -> update_by(acceleration, step));
    retval.push_back(position -> update_by(&retval[0], step));
    return retval;
}

template <>
vec2 rk4(double m1, double m2, vec2 v1, vec2 v2, vec2 p1, vec2 p2, Simulation<vec2>* _s){
    vec2 k1, k2, k3, k4;
    std::vector<vec2> dummy;
    k1 = _s -> get_force(m1, m2, v1, v2, p1, p2, _s) * _s -> get_dt() ;
    dummy = get_new_pos(&p1, &v1, &k1, _s -> get_dt() /2);
    k2 = _s -> get_force(m1, m2, v1 + dummy[0]/2.f, v2, p1 + dummy[1]/2.f, p2, _s) * _s -> get_dt() ;
    dummy = get_new_pos(&p1, &v1, &k2, _s -> get_dt() /2);
    k3 = _s -> get_force(m1, m2, v1 + dummy[0]/2.f, v2, p1 + dummy[1]/2.f, p2, _s) * _s -> get_dt() ;
    dummy = get_new_pos(&p1, &v1, &k3, _s -> get_dt() );
    k4 = _s -> get_force(m1, m2, v1 + dummy[0]/2.f, v2, p1 + dummy[1]/2.f, p2, _s) * _s -> get_dt() ;

    return (k1 + k2*2.f + k3*2.f + k4)/6.f;
}

//===---------------------------------------------------------===//
// 3d Force calculation methods                                  //
//===---------------------------------------------------------===//

template <>
void update_euler(Body<vec3>* b, Simulation<vec3>* _s){
    b -> velocity += b -> acceleration* _s -> get_dt() ;
    b -> position += b -> velocity * _s -> get_dt() ;
    b -> acceleration.reset();
}

template <>
void update_leapfrog(Body<vec3>* b, Simulation<vec3>* _s){
   b -> position += b -> velocity * _s -> get_dt()  + b -> acceleration * _s -> get_dt()  * _s -> get_dt()  * .5;
   b -> velocity += (b -> acceleration + b -> prev_acc)* _s -> get_dt()  * .5;
   b -> acceleration.reset();
}

template <>
void update_symplectic4(Body<vec3>* body, Simulation<vec3>* _s){
    body -> position += body -> velocity * c1 * _s -> get_dt() ;
    
    _s -> update_pair(body);
    body -> velocity += body -> acceleration * d1 *_s -> get_dt();
    body -> position += body -> velocity * c2 *_s -> get_dt();
    
    _s -> update_pair(body);
    body -> velocity += body -> acceleration * d2 *_s -> get_dt();
    body -> position += body -> velocity * c2 *_s -> get_dt();
    
    _s -> update_pair(body);
    body -> velocity += body -> acceleration * d1 *_s -> get_dt();
    body -> position += body -> velocity * c1 * _s -> get_dt() ;
}

//===---------------------------------------------------------===//
// Force calculation methods                                     //
//===---------------------------------------------------------===//

template <>
vec3 newtonian(double m1, double m2, vec3 v1, vec3 v2, vec3 p1, vec3 p2, Simulation<vec3>* _s){
    vec3 d = p2 - p1;
    return d.normalize() *_s -> get_G()* m1*m1/ (d.sqr_magn()+ 1);
}

template <>
vec3 pn2(double m1, double m2, vec3 v1, vec3 v2, vec3 p1, vec3 p2, Simulation<vec3>* _s){ 
    vec3 n = p2 - p1;
    double r = n.magnitude() * _s -> get_scale()+1;

    double term1 = _s -> get_G()*m2/(_s -> get_c_sqr()*r*r);
    double m_ratio = m1/m2;
    vec3 v = v1 - v2;

    double correction = std::exp(-r*r/_s -> get_takeover())*_s -> get_e_sqr();
    vec3 a_newton = n * _s -> get_G() *m1*m2/(r*r);
    vec3 pn1 = ( n*(4 + 2*m_ratio)*term1*_s -> get_c_sqr() - v*n*(1+ 3*m_ratio) + v*n*v*4 - v*n*n*3/2) * term1 ;
    vec3 pn25 = (v*n)*v*-8/5 * _s -> get_G()* _s -> get_G()*m1*m2*(m1+m2)/(_s -> get_c_sqr() * _s -> get_c_sqr()* _s -> get_c()*r*r*r);

    return a_newton + pn1 + pn25 + n*correction;
}

template <>
std::vector<vec3> get_new_pos(vec3* position, vec3* velocity, vec3* acceleration, double step){
    std::vector<vec3> retval;
    retval.push_back(velocity -> update_by(acceleration, step));
    retval.push_back(position -> update_by(&retval[0], step));
    return retval;
}

template <>
vec3 rk4(double m1, double m2, vec3 v1, vec3 v2, vec3 p1, vec3 p2, Simulation<vec3>* _s){
    vec3 k1, k2, k3, k4;
    std::vector<vec3> dummy;
    k1 = _s -> get_force(m1, m2, v1, v2, p1, p2, _s) * _s -> get_dt() ;
    dummy = get_new_pos(&p1, &v1, &k1, _s -> get_dt() /2);
    k2 = _s -> get_force(m1, m2, v1 + dummy[0]/2.f, v2, p1 + dummy[1]/2.f, p2, _s) * _s -> get_dt() ;
    dummy = get_new_pos(&p1, &v1, &k2, _s -> get_dt() /2);
    k3 = _s -> get_force(m1, m2, v1 + dummy[0]/2.f, v2, p1 + dummy[1]/2.f, p2, _s) * _s -> get_dt() ;
    dummy = get_new_pos(&p1, &v1, &k3, _s -> get_dt() );
    k4 = _s -> get_force(m1, m2, v1 + dummy[0]/2.f, v2, p1 + dummy[1]/2.f, p2, _s) * _s -> get_dt() ;

    return (k1 + k2*2.f + k3*2.f + k4)/6.f;
}
}
