#pragma once

#include <cmath>

#include "Aster/physics/tool-chain.h"
#include "Aster/physics/vectors.h"

#include "Aster/simulations/sim_obj.h"
#include "Aster/simulations/basic.h"

namespace Aster{


double c1 = 1 / (2 * (2 - std::pow(2.0, 1.0/3.0)));
double c2 = (1 - pow(2.0, 1.0/3.0)) * c1;
double d1 = 1 / (2 - std::pow(2.0, 1.0/3.0));
double d2 = -std::pow(2.0, 1.0/3.0) *d1; 
//===---------------------------------------------------------===//
// Update methods                                                //
//===---------------------------------------------------------===//


extern const double PI;

template <typename T> 
inline double get_A1(double eta, T v, double r_dot, double m, double r){
    // -(1 + 3η)v² + 1.5ηv² + 2(2 + η)m/r
    return - (1 + 3*eta) * v.sqr_magn() + 3/2 * eta * r_dot * r_dot + 2*(2 + eta) * m / r; 
}

template <typename T> 
inline double get_B1(double eta){
    // 2(2 - η)
    return 2*(2 - eta); 
}

template <typename T> 
inline double get_A2(double eta, T v, double r_dot, double m, double r){
    // -η(3 - 4η)v^4 + .5η(13 - 4η)v²m/r + 3/2η(3 - 4η)v²r_dot² + (2 + 25η + 2η²)r_dot²m/r - 15/8η(1 - 3η)r_dot^4 - 3/4 (12 + 29η)(m/r)²
    double retval = - eta * (3 - 4 * eta) * std::pow(v.sqr_magn(), 2);
    retval += 1.0/2.0 * eta * (13 - 4*eta) * v.sqr_magn() * m / r;
    retval += 3.0/2.0 * eta * (3 - 4*eta) * v.sqr_magn() * r_dot * r_dot;
    retval += (2 + 25*eta + 2 *eta *eta) * r_dot * r_dot * m /r;
    retval += - 15.0/8.0 * eta * (1 - 3*eta) * std::pow(r_dot, 4) - 3.0/4.0 * (12 + 29 * eta) * (m /r) * (m/r);

    return retval;
}

template <typename T> 
inline double get_B2(double eta, T v, double r_dot, double m, double r){
    double retval = 1/2 * eta * (15 + 4 *eta) * v.sqr_magn();
    retval += -3.0/2.0 * eta * (3+ 2*eta) * r_dot *r_dot;
    retval += -1.0/2 * (4 + 41*eta + 8*eta*eta) * m /r;

    return retval;
}

template <typename T> 
inline double get_A25(double eta, T v, double r_dot, double m, double r){
    return 3 * v.sqr_magn() + 17.0/3.0 * m /r;
}

template <typename T> 
inline double get_B25(double eta, T v, double r_dot, double m, double r){
    return v.sqr_magn() + 3 * m /r;
}

template <typename T> 
inline T pn25(double m1, double m2, T v1, T v2, T p1, T p2, Simulation<T>* _s){
    T x = p1 - p2;
    double r = x.magnitude() + 1;
    double m = m1 + m2;
    double eta = m1*m2 / (m*m);
    
    T v = v1 - v2;
    T n = x / r;
    double r_dot = v * n;

    double a_components = get_A1<T>(eta, v, r_dot, m, r) + get_A2<T>(eta, v, r_dot, m, r);
    double b_components = get_B1<T>(eta                ) + get_B2<T>(eta, v, r_dot, m, r);

    double half_a = get_A25<T>(eta, v, r_dot, m, r);
    double half_b = get_B25<T>(eta, v, r_dot, m, r);

    T acc = - n* m / (r*r); 
    acc +=  ((n * a_components + v *r_dot * b_components) * m / (r*r)) / std::pow(_s -> get_c(), 3); 
    acc += ( 8.0/5.0 * eta * (m*m) / (r*r*r) * (n * r_dot * half_a - v * half_b)) / (std::pow(_s -> get_c(), 5));

    return acc * _s -> get_G() * m1;
}

template <typename T> 
inline T pn2(double m1, double m2, T v1, T v2, T p1, T p2, Simulation<T>* _s){
    T x = p1 - p2;
    double r = x.magnitude() + 1;
    double m = m1 + m2;
    double eta = m1*m2 / (m*m);
    
    T v = v1 - v2;
    T n = x / r;
    double r_dot = v * n;

    double a_components = get_A1<T>(eta, v, r_dot, m, r) + get_A2<T>(eta, v, r_dot, m, r);
    double b_components = get_B1<T>(eta                ) + get_B2<T>(eta, v, r_dot, m, r);


    T acc = - n* m / (r*r); 
    acc +=  ((n * a_components + v *r_dot * b_components) * m / (r*r)) / std::pow(_s -> get_c(), 3);  
    return acc* _s -> get_G() * m1;
}

template <typename T> 
inline T pn1(double m1, double m2, T v1, T v2, T p1, T p2, Simulation<T>* _s){
    T x = p1 - p2;
    double r = x.magnitude() + 1;
    double m = m1 + m2;
    double eta = m1*m2 / (m*m);
    
    T v = v1 - v2;
    T n = x / r;
    double r_dot = v * n;

    double a_components = get_A1<T>(eta, v, r_dot, m, r);
    double b_components = get_B1<T>(eta                );


    T acc = - n * m / (r*r); 
    acc +=  ((n * a_components + v *r_dot * b_components) * m / (r*r)) / std::pow(_s -> get_c(), 3);  
    return acc* _s -> get_G() * m1;
}



template <typename T>
double get_eccentricity(Simulation<T>* _s, Body<T>* body, double relv_sq, double w_squared,  double radius, double mass2){
    double mu = _s -> get_G() * (body -> mass + mass2);
    double reduced_mass  = body -> mass * mass2 / (body -> mass + mass2 ); 
    double orbital_energy = relv_sq / 2 - mu / radius;

    return 1 + 2 * orbital_energy * w_squared / (reduced_mass * mu * mu + _s -> get_e_sqr());
}



template <typename T>
void compute_rad_pressure(Simulation<T>* _s, Body<T>* body, T pos, double temp){
    T vec = (pos - body -> position);
    double r_sq = vec.sqr_magn();
    double force = _s -> get_boltzmann() * std::pow(body -> temp, 4) / (4 * PI * r_sq * _s -> get_c());

    body -> acceleration += - vec.normalize() * force / body -> mass;
}

template <>
void compute_rad_pressure(Simulation<vec2>* _s, Body<vec2>* body, vec2 pos, double temp){
    vec2 vec = (pos - body -> position);
    double r_sq = vec.sqr_magn();
    double force = _s -> get_boltzmann() * std::pow(body -> temp, 4) / (4 * PI * r_sq * _s -> get_c());

    body -> acceleration += - vec.normalize() * force / body -> mass;
}

template <>
void get_new_temp<vec3>(Simulation<vec3>* _s, Body<vec3>* body, vec3 pos, vec3 vel, double temp, double mass){ 
    /*
    * we are trying to compute the delta T so we use the formula 
    *  21 * G * M * m * n * e²
    *  ------------------------ = E_t
    *         2 * r^6

    * from here it is trival, we can just divide by the mass and a coeff and multiply 
    * byt dt to get the delta T 
    */
    // sketchy edge case...
    if (body -> acceleration.sqr_magn()) return;

    vec3 rel_v = body -> velocity - vel;
    double relv_sq = rel_v.sqr_magn();
    double r = (body -> position - pos).magnitude();
    
    /*
    * omega is the angular velocity of the body
    * thus it can be calculated by getting the perpendicular
    * vector to the projection of the velocity onto the acceleration
    * assuming the acceleration comes purely from the mutual
    * attraction between the bodies
    * ergo 
    * W  = V x a / a x a * vec(a)
    */
    double w_squared = (body -> acceleration * (rel_v * body -> acceleration) / (body -> acceleration.sqr_magn())).sqr_magn();

    double eccentricity = get_eccentricity<vec3>(_s, body, relv_sq, w_squared, r, mass);

    double delta_temperature = 1;//21.0 / 2.0 * _s -> get_G() * mass * body -> mass * eccentricity * eccentricity * w_squared / (r*r*r*r*r*r);
    delta_temperature *= _s -> get_dt(); // we want it to be per unit of time
    delta_temperature /= body -> mass * _s -> get_heat_capacity(); // then we get how many K the body got in dt 

    // boltzmann
    delta_temperature -= _s -> get_boltzmann() * std::pow(body -> temp, 4) / ( _s -> get_heat_capacity() * body -> mass);

    body -> temp += delta_temperature; // done! 
}

template <>
void get_new_temp<vec2>(Simulation<vec2>* _s, Body<vec2>* body, vec2 pos, vec2 vel, double temp, double mass){
    /*
    * we are trying to compute the delta T so we use the formula 
    *  21 * G * M * m * n * e²
    *  ------------------------ = E_t
    *         2 * r^6

    * from here it is trival, we can just divide by the mass and a coeff and multiply 
    * byt dt to get the delta T 
    */

    vec2 rel_v = body -> velocity - vel;
    double relv_sq = rel_v.sqr_magn();
    double r = (body -> position - pos).magnitude();
    /*
    * omega^2 = |V|^2 - (V x A / |A|) ^2 
    * omega being the hypoteetase of a right triangle
    * of sidelenghts equal to the velocity and the
    * projection of the velocity onto the acceleration
    */
    double w_squared = (body -> acceleration * (rel_v * body -> acceleration) / body -> acceleration.sqr_magn()).sqr_magn();

    double eccentricity = get_eccentricity<vec2>(_s, body, relv_sq, w_squared, r, mass);

    double delta_temperature = 21.0 / 2.0 * _s -> get_G() * mass * body -> mass * eccentricity * eccentricity * w_squared / (r*r*r*r*r*r);
    delta_temperature *= _s -> get_dt(); // we want it to be per unit of time
    delta_temperature /= body -> mass * _s -> get_heat_capacity(); // then we get how many K the body got in dt 

    // boltzmann
    delta_temperature -= _s -> get_boltzmann() * std::pow(body -> temp, 4) / ( _s -> get_heat_capacity() * body -> mass);

    body -> temp +=  delta_temperature; // done! 
}



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

template <typename T>
T newtonian(double m1, double m2, T v1, T v2, T p1, T p2, Simulation<T>* _s){
    T d = p2 - p1;
    return d.normalize() *_s -> get_G()* m1*m1/ (d.sqr_magn()+ 1);
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
    //vec3 pn1 = ( n*(4 + 2*m_ratio)*term1*_s -> get_c_sqr() - v*n*(1+ 3*m_ratio) + v*n*v*4 - v*n*n*3/2) * term1 ;
    //vec3 pn25 = (v*n)*v*-8/5 * _s -> get_G()* _s -> get_G()*m1*m2*(m1+m2)/(_s -> get_c_sqr() * _s -> get_c_sqr()* _s -> get_c()*r*r*r);

    return a_newton;// + pn1 + pn25 + n*correction;
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

template <>
func_ptr<vec2> get_update_func(update_type type){
    switch (type){
    case EULER:
        return update_euler<vec2>;
    case LEAPFROG:
        return update_leapfrog<vec2>;
    case SYMPLECTIC4:
        return update_symplectic4<vec2>;
    default:
        return update_euler<vec2>;
    }
}

template <>
force_func<vec2> get_force_func(force_type type){
    switch (type){
    case NEWTON:
        return newtonian<vec2>;
    case PN1:
        return pn1<vec2>;
    case PN2:
        return pn2<vec2>;
    case PN25:
        return pn25<vec2>;
    default:
        return newtonian<vec2>;
    }   
}

template <>
func_ptr<vec3> get_update_func(update_type type){
    switch (type){
    case EULER:
        return update_euler<vec3>;
    case LEAPFROG:
        return update_leapfrog<vec3>;
    case SYMPLECTIC4:
        return update_symplectic4<vec3>;
    default:
        return update_euler<vec3>;
    }
}

template <>
force_func<vec3> get_force_func(force_type type){
    switch (type){
    case NEWTON:
        return newtonian<vec3>;
    case PN1:
        return pn1<vec3>;
    case PN2:
        return pn2<vec3>;
    case PN25:
        return pn25<vec3>;
    default:
        return newtonian<vec3>;
    }   
}
}
