#pragma once

#include <cmath>

#include "Aster/physics/tool-chain.h"
#include "Aster/physics/vectors.h"

#include "Aster/simulations/sim_obj.h"
#include "Aster/simulations/basic.h"

#include "Aster/impl/config.h"

namespace Aster{ 
//===---------------------------------------------------------===//
// Update methods                                                //
//===---------------------------------------------------------===//

// !!! POSITIVE FORCE = ATTRACTION

extern const double PI;

template <typename T, typename F>
FORCE_INLINE void for_each_body(Simulation<T>* _s, F func) {
    assert(_s->has_loaded_yet());
    
    const size_t n_bodies = _s->bodies.size();
    const size_t n_threads = std::min(n_bodies, static_cast<size_t>(_s->get_cores()));
    
    if (n_bodies == 0) return;
    
    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    
    // Calcola la dimensione dei chunk
    const size_t chunk_size = (n_bodies + n_threads - 1) / n_threads;
    
    // Crea i thread
    for (size_t i = 0; i < n_threads; ++i) {
        const size_t start = i * chunk_size;
        const size_t end = std::min(start + chunk_size, n_bodies);
        
        if (start >= end) break;
        
        threads.emplace_back([start, end, &func, _s]() {
            for (size_t b = start; b < end; ++b) {
                func(_s->bodies[b]);
            }
        });
    }
    
    // Attendi il completamento
    for (auto& t : threads) {
        t.join();
    }
}

template <typename T> 
FORCE_INLINE double get_A1(double eta, T v, double r_dot, double m, double r){
    // -(1 + 3η)v² + 1.5ηv² + 2(2 + η)m/r
    return - (1 + 3*eta) * v.sqr_magn() + 3/2 * eta * r_dot * r_dot + 2*(2 + eta) * m / r; 
}

template <typename T> 
FORCE_INLINE double get_B1(double eta){
    // 2(2 - η)
    return 2*(2 - eta); 
}

template <typename T> 
FORCE_INLINE double get_A2(double eta, T v, double r_dot, double m, double r){
    // -η(3 - 4η)v^4 + .5η(13 - 4η)v²m/r + 3/2η(3 - 4η)v²r_dot² + (2 + 25η + 2η²)r_dot²m/r - 15/8η(1 - 3η)r_dot^4 - 3/4 (12 + 29η)(m/r)²
    double retval = - eta * (3 - 4 * eta) * std::pow(v.sqr_magn(), 2);
    retval += 1.0/2.0 * eta * (13 - 4*eta) * v.sqr_magn() * m / r;
    retval += 3.0/2.0 * eta * (3 - 4*eta) * v.sqr_magn() * r_dot * r_dot;
    retval += (2 + 25*eta + 2 *eta *eta) * r_dot * r_dot * m /r;
    retval += - 15.0/8.0 * eta * (1 - 3*eta) * std::pow(r_dot, 4) - 3.0/4.0 * (12 + 29 * eta) * (m /r) * (m/r);

    return retval;
}

template <typename T> 
FORCE_INLINE double get_B2(double eta, T v, double r_dot, double m, double r){
    double retval = 1/2 * eta * (15 + 4 *eta) * v.sqr_magn();
    retval += -3.0/2.0 * eta * (3+ 2*eta) * r_dot *r_dot;
    retval += -1.0/2 * (4 + 41*eta + 8*eta*eta) * m /r;

    return retval;
}

template <typename T> 
FORCE_INLINE double get_A25(double eta, T v, double r_dot, double m, double r){
    return 3 * v.sqr_magn() + 17.0/3.0 * m /r;
}

template <typename T> 
FORCE_INLINE double get_B25(double eta, T v, double r_dot, double m, double r){
    return v.sqr_magn() + 3 * m /r;
}

template <typename T> 
inline  T pn25(double m1, double m2, T v1, T v2, T p1, T p2, Simulation<T>* _s){
    T x = (p2 - p1 + T(_s -> get_e_sqr()));
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

    T acc = n * m1*m2 / (r*r); 
    acc +=  -((n * a_components + v *r_dot * b_components) * m / (r*r)) / std::pow(_s -> get_c(), 3); 
    acc += -( 8.0/5.0 * eta * (m*m) / (r*r*r) * (n * r_dot * half_a - v * half_b)) / (std::pow(_s -> get_c(), 5));

    return acc * _s -> get_G();
}

template <typename T> 
inline T pn2(double m1, double m2, T v1, T v2, T p1, T p2, Simulation<T>* _s){
    T x = (p2 - p1 + T(_s -> get_e_sqr()));
    double r = x.magnitude();
    double m = m1 + m2;
    double eta = m1*m2 / (m*m);
    
    T v = v1 - v2;
    T n = x / r;
    double r_dot = v * n;

    double a_components = get_A1<T>(eta, v, r_dot, m, r) + get_A2<T>(eta, v, r_dot, m, r);
    double b_components = get_B1<T>(eta                ) + get_B2<T>(eta, v, r_dot, m, r);


    T acc = n * m1*m2 / (r*r); 
    acc +=  -((n * a_components + v *r_dot * b_components) * m / (r*r)) / std::pow(_s -> get_c(), 3);  
    return acc* _s -> get_G();
}

template <typename T> 
inline T pn1(double m1, double m2, T v1, T v2, T p1, T p2, Simulation<T>* _s){
    T x = (p2 - p1 + T(_s -> get_e_sqr()));
    double r = x.magnitude();
    double m = m1 + m2;
    double eta = m1*m2 / (m*m);
    
    T v = v1 - v2;
    T n = x / r;
    double r_dot = v * n;

    double a_components = get_A1<T>(eta, v, r_dot, m, r);
    double b_components = get_B1<T>(eta                );


    T acc = n * m1*m2 / (r*r); 
    acc +=  -((n * a_components + v *r_dot * b_components) * m / (r*r)) / std::pow(_s -> get_c(), 3);  
    return acc* _s -> get_G();
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
void update_euler(Simulation<vec2>* _s){
    for (auto& body : _s -> bodies){
        body.velocity += body.acceleration* _s -> get_dt() ;
        body.position += body.velocity * _s -> get_dt() ;
    }
}

template <>
void update_leapfrog(Simulation<vec2>* _s){
    for (auto& body : _s -> bodies){
        body.position += body.velocity * _s -> get_dt()  + body.acceleration * _s -> get_dt()  * _s -> get_dt()  * .5;
        body.velocity += (body.acceleration + body.prev_acc)* _s -> get_dt()  * .5;
    }
}

template <>
void update_symplectic4(Simulation<vec2>* _s){
    constexpr double c1 = 0.6756035959798289;
    constexpr double c2 = -0.17560359597982883;
    constexpr double d1 = 1.3512071919596578;
    constexpr double d2 = -1.7024143839193153;

    static vec2 temp_v;

    for (auto& body : _s -> bodies){
        // first step 
        temp_v = body.velocity + d1*body.acceleration * _s -> get_dt();

        body.position += c1* temp_v * _s -> get_dt();
        body.velocity = temp_v;

        // second
        _s -> update_forces();
        temp_v = body.velocity + d2*body.acceleration * _s -> get_dt();

        body.position += c2* temp_v * _s -> get_dt();
        body.velocity = temp_v;

        //third
        _s -> update_forces();
        temp_v = body.velocity + d1*body.acceleration * _s -> get_dt();

        body.position += c2* temp_v * _s -> get_dt();
        body.velocity = temp_v;

        //last
        _s -> update_forces();
        body.position += c1* body.velocity * _s -> get_dt();
    }
}

template <typename T> 
void adaptive_euler(Simulation<T>* _s){
    for (auto& body : _s -> bodies){
        double magn = body.acceleration.magnitude();
        T r = body.acceleration / magn;

        double E_i = body.velocity.sqr_magn() * .5 *body.mass - std::sqrt(_s -> get_G() * (_s -> get_total_mass() - body.mass) / magn);

        Body<T> rollback = body;
        update_euler<T>(_s);

        double E_f = body.velocity.sqr_magn() * .5 *body.mass - std::sqrt(_s -> get_G() * (_s -> get_total_mass() - body.mass) / body.acceleration.magnitude());
        double iters = std::log10(std::abs((E_i - E_f) / E_i));

        if (iters < _s -> get_adaptive_coeff())
            return;

        iters = _s -> get_adaptive_coeff() - iters;

        body.position = rollback.position;
        body.acceleration = rollback.acceleration;
        body.temp = rollback.temp;
        body.velocity = rollback.velocity;

        _s -> set_dt(_s -> get_dt() / iters);
        for (double i = 0; i < iters; ++i){
            update_euler<T>(_s);
        }

        _s -> set_dt(_s -> get_dt() * iters);
    }
}
//===---------------------------------------------------------===//
// Force calculation methods                                     //
//===---------------------------------------------------------===//

template <typename T>
T newtonian(double m1, double m2, T v1, T v2, T p1, T p2, Simulation<T>* _s){
    T n = (p2 - p1);
    double r = n.magnitude();

    return n.normalize() * _s -> get_G() *m1*m2/(r*r);
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
// 3d Force update methods                                       //
//===---------------------------------------------------------===//

template <>
void update_euler(Simulation<vec3>* _s){
    for_each_body(_s, [_s](Body<vec3>& body){
        body.velocity += body.acceleration* _s -> get_dt() ;
        body.position += body.velocity * _s -> get_dt() ;
    });
}

template <>
void update_leapfrog(Simulation<vec3>* _s){
    for_each_body(_s, [_s](Body<vec3>& body){
        body.position += body.velocity * _s -> get_dt()  + body.acceleration * _s -> get_dt()  * _s -> get_dt()  * .5;
        body.velocity += (body.acceleration + body.prev_acc)* _s -> get_dt()  * .5;
    });
}

template <>
void update_symplectic4(Simulation<vec3>* _s){
    constexpr double c1 = 0.675603595979828885909057589742587879;
    constexpr double c2 = -0.175603595979828858153481974113674369;
    constexpr double d1 = 1.35120719195965777181811517948517576;
    constexpr double d2 = -1.70241438391931532159162543393904343;
    
    // first step 
    for_each_body(_s, [d1, c1, _s](Body<vec3>& body){
        vec3 temp_v = body.velocity + d1*body.acceleration * _s -> get_dt();

        body.position += c1* temp_v * _s -> get_dt();
        body.velocity = temp_v;
    });

    // second
    _s -> update_forces();

    for_each_body(_s, [c2, _s](Body<vec3>& body){
        vec3 temp_v = body.velocity + d2*body.acceleration * _s -> get_dt();

        body.position += c2* temp_v * _s -> get_dt();
        body.velocity = temp_v;
    });
    
    //third
    _s -> update_forces();

    for_each_body(_s, [d1, c2, _s](Body<vec3>& body){
        vec3 temp_v = body.velocity + d1*body.acceleration * _s -> get_dt();

        body.position += c2* temp_v * _s -> get_dt();
        body.velocity = temp_v;
    });

    //last
    _s -> update_forces();

    for_each_body(_s, [c1, _s](Body<vec3>& body){
        body.position += c1 * body.velocity * _s -> get_dt();
    });
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
    case ADE:
        return adaptive_euler<vec2>;
    case SABA1:
        return update_SABA1<vec2>;
    case SABA2:
        return update_SABA2<vec2>;
    case SABA3:
        return update_SABA3<vec2>;
    case SABA4:
        return update_SABA4<vec2>;
    case SABA5:
        return update_SABA5<vec2>;
    case SABA6:
        return update_SABA6<vec2>;
    case SABA7:
        return update_SABA7<vec2>;
    case SABA8:
        return update_SABA8<vec2>;
    case SABA9:
        return update_SABA9<vec2>;
    case SABA10:
        return update_SABA10<vec2>;
    default:
        warn_if(true, "update method not found");
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
        warn_if(true, "force method not found");
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
    case ADE:
        return adaptive_euler<vec3>;
    case SABA1:
        return update_SABA1<vec3>;
    case SABA2:
        return update_SABA2<vec3>;
    case SABA3:
        return update_SABA3<vec3>;
    case SABA4:
        return update_SABA4<vec3>;
    case SABA5:
        return update_SABA5<vec3>;
    case SABA6:
        return update_SABA6<vec3>;
    case SABA7:
        return update_SABA7<vec3>;
    case SABA8:
        return update_SABA8<vec3>;
    case SABA9:
        return update_SABA9<vec3>;
    case SABA10:
        return update_SABA10<vec3>;
    default:
        warn_if(true, "update method not found");
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
        warn_if(true, "force method not found");
        return newtonian<vec3>;
    }   
}
}