#include <cmath>

#include <tbb/parallel_for.h>
#include "Aster/physics/vectors.h"

#include "Aster/simulations/sim_obj.h"
#include "Aster/physics/tool-chain.h"
#include "Aster/simulations/basic.h"
#include "Aster/building-api/logging.h"

namespace Aster{ 
using namespace tbb;


//===---------------------------------------------------------===//
// Update methods                                                //
//===---------------------------------------------------------===//

//! POSITIVE FORCE = ATTRACTION

extern const REAL PI;

 
FORCE_INLINE REAL get_A1(REAL eta, vec3 v, REAL r_dot, REAL m, REAL r){
    // -(1 + 3η)v² + 1.5ηv² + 2(2 + η)m/r
    return - (1 + 3*eta) * v.sqr_magn() + 3/2 * eta * r_dot * r_dot + 2*(2 + eta) * m / r; 
}

 
FORCE_INLINE REAL get_B1(REAL eta){
    // 2(2 - η)
    return 2*(2 - eta); 
}

 
FORCE_INLINE REAL get_A2(REAL eta, vec3 v, REAL r_dot, REAL m, REAL r){
    // -η(3 - 4η)v^4 + .5η(13 - 4η)v²m/r + 3/2η(3 - 4η)v²r_dot² + (2 + 25η + 2η²)r_dot²m/r - 15/8η(1 - 3η)r_dot^4 - 3/4 (12 + 29η)(m/r)²
    REAL retval = - eta * (3 - 4 * eta) * std::pow(v.sqr_magn(), 2);
    retval += 1.0/2.0 * eta * (13 - 4*eta) * v.sqr_magn() * m / r;
    retval += 3.0/2.0 * eta * (3 - 4*eta) * v.sqr_magn() * r_dot * r_dot;
    retval += (2 + 25*eta + 2 *eta *eta) * r_dot * r_dot * m /r;
    retval += - 15.0/8.0 * eta * (1 - 3*eta) * std::pow(r_dot, 4) - 3.0/4.0 * (12 + 29 * eta) * (m /r) * (m/r);

    return retval;
}

 
FORCE_INLINE REAL get_B2(REAL eta, vec3 v, REAL r_dot, REAL m, REAL r){
    REAL retval = 1/2 * eta * (15 + 4 *eta) * v.sqr_magn();
    retval += -3.0/2.0 * eta * (3+ 2*eta) * r_dot *r_dot;
    retval += -1.0/2 * (4 + 41*eta + 8*eta*eta) * m /r;

    return retval;
}

 
FORCE_INLINE REAL get_A25(REAL eta, vec3 v, REAL r_dot, REAL m, REAL r){
    return 3 * v.sqr_magn() + 17.0/3.0 * m /r;
}

 
FORCE_INLINE REAL get_B25(REAL eta, vec3 v, REAL r_dot, REAL m, REAL r){
    return v.sqr_magn() + 3 * m /r;
}

 
vec3 pn25(REAL m1, REAL m2, vec3 v1, vec3 v2, vec3 p1, vec3 p2, Simulation* _s){
    vec3 x = (p2 - p1 + vec3(_s -> get_e_sqr()));
    REAL r = x.magnitude() + 1;
    REAL m = m1 + m2;
    REAL eta = m1*m2 / (m*m);
    
    vec3 v = v1 - v2;
    vec3 n = x / r;
    REAL r_dot = v * n;

    REAL a_components = get_A1(eta, v, r_dot, m, r) + get_A2(eta, v, r_dot, m, r);
    REAL b_components = get_B1(eta                ) + get_B2(eta, v, r_dot, m, r);

    REAL half_a = get_A25(eta, v, r_dot, m, r);
    REAL half_b = get_B25(eta, v, r_dot, m, r);

    vec3 acc = n * m1*m2 / (r*r+ _s->softening); 
    acc +=  -((n * a_components + v *r_dot * b_components) * m / (r*r+ _s->softening)) / std::pow(_s -> get_c(), 3); 
    acc += -( 8.0/5.0 * eta * (m*m) / (r*r*r  + _s->softening) * (n * r_dot * half_a - v * half_b)) / (std::pow(_s -> get_c(), 5));

    return acc * _s -> get_G();
}

 
vec3 pn2(REAL m1, REAL m2, vec3 v1, vec3 v2, vec3 p1, vec3 p2, Simulation* _s){
    vec3 x = (p2 - p1 + vec3(_s -> get_e_sqr()));
    REAL r = x.magnitude();
    REAL m = m1 + m2;
    REAL eta = m1*m2 / (m*m);
    
    vec3 v = v1 - v2;
    vec3 n = x / r;
    REAL r_dot = v * n;

    REAL a_components = get_A1(eta, v, r_dot, m, r) + get_A2(eta, v, r_dot, m, r);
    REAL b_components = get_B1(eta                ) + get_B2(eta, v, r_dot, m, r);


    vec3 acc = n * m1*m2 / (r*r + _s->softening); 
    acc +=  -((n * a_components + v *r_dot * b_components) * m / (r*r + _s->softening)) / std::pow(_s -> get_c(), 3);  
    return acc* _s -> get_G();
}

 
vec3 pn1(REAL m1, REAL m2, vec3 v1, vec3 v2, vec3 p1, vec3 p2, Simulation* _s){
    vec3 x = (p2 - p1 + vec3(_s -> get_e_sqr()));
    REAL r = x.magnitude();
    REAL m = m1 + m2;
    REAL eta = m1*m2 / (m*m);
    
    vec3 v = v1 - v2;
    vec3 n = x / r;
    REAL r_dot = v * n;

    REAL a_components = get_A1(eta, v, r_dot, m, r);
    REAL b_components = get_B1(eta                );


    vec3 acc = n * m1*m2 / (r*r+ _s->softening); 
    acc +=  -((n * a_components + v *r_dot * b_components) * m / (r*r+ _s->softening)) / std::pow(_s -> get_c(), 3);  
    return acc* _s -> get_G();
}



REAL get_eccentricity(Simulation* _s, size_t body, REAL relv_sq, REAL w_squared,  REAL radius, REAL mass2){
    REAL mu = _s -> get_G() * (_s -> bodies.get_mass_of(body) + mass2);
    REAL reduced_mass  = _s -> bodies.get_mass_of(body) * mass2 / (_s -> bodies.get_mass_of(body) + mass2 ); 
    REAL orbital_energy = relv_sq / 2 - mu / radius;

    return 1 + 2 * orbital_energy * w_squared / (reduced_mass * mu * mu + _s -> get_e_sqr());
}


inline void compute_rad_pressure(Simulation* _s, size_t body, vec3 pos, REAL temp){
    vec3 vec = (pos - _s -> bodies.get_position_of(body));
    REAL r_sq = vec.sqr_magn();
    REAL force = _s -> get_boltzmann() * std::pow(_s -> bodies.get_temp_of(body), 4) / (4 * PI * r_sq * _s -> get_c());

    _s -> bodies.get_acc_of(body) += - vec.normalize() * force / _s -> bodies.get_mass_of(body);
}


void update_euler(Simulation* _s){
    _s -> solver -> compute_forces();

    size_t N = _s -> bodies.positions.size();
    parallel_for(size_t(0), N, [_s](size_t i){
        _s -> bodies.velocities[i] += _s -> bodies.accs[i] * _s -> get_dt();
        _s -> bodies.positions[i] += _s -> bodies.velocities[i] * _s -> get_dt();
        _s -> bodies.accs[i].reset();
    });
}

void update_leapfrog(Simulation* _s){
    _s -> solver -> compute_forces();

    size_t N = _s -> bodies.positions.size();
    REAL dt = _s -> get_dt();

    parallel_for(size_t(0), N, [&, _s, dt](size_t i){
        _s -> bodies.positions[i] += _s -> bodies.velocities[i] * dt + .5 *_s -> bodies.accs[i] * dt * dt; 
        _s -> bodies.velocities[i] += .5 * (_s -> bodies.accs[i]) * dt;
    });

    _s -> solver -> compute_forces();

    parallel_for(size_t(0), N, [&, _s, dt](size_t i){
        _s -> bodies.velocities[i] += .5 * (_s -> bodies.accs[i]) * dt; 
    });
}

inline void update_kepler(vec3& r0, vec3& v0, REAL mu, REAL dt) {
    constexpr REAL tol = 1e-16; 
    constexpr int maxiter = 500;
    constexpr REAL z_series_thresh = 1e-8;

    auto stumpff_C = [&](REAL z) -> REAL {
        if (std::abs(z) < z_series_thresh) {
            REAL z2 = z*z;
            REAL z3 = z2*z;
            return 0.5 - z/24.0 + z2/720.0 - z3/40320.0;
        } else if (z > 0.0) {
            REAL s = std::sqrt(z);
            return (1.0 - std::cos(s)) / z;
        } else {
            REAL s = std::sqrt(-z);
            return (std::cosh(s) - 1.0) / (-z);
        }
    };

    auto stumpff_S = [&](REAL z) -> REAL {
        if (std::abs(z) < z_series_thresh) {
            REAL z2 = z*z;
            REAL z3 = z2*z;
            return (1.0/6.0) - z/120.0 + z2/5040.0 - z3/362880.0;
        } else if (z > 0.0) {
            REAL s = std::sqrt(z);
            return (s - std::sin(s)) / (s*s*s);
        } else {
            REAL s = std::sqrt(-z);
            return (std::sinh(s) - s) / (s*s*s);
        }
    };
    REAL r0n = r0.magnitude();
    REAL v0n = v0.magnitude();
    REAL vr0 = (r0.x * v0.x + r0.y * v0.y + r0.z * v0.z) / r0n;

    REAL energy = v0n*v0n/2.0 - mu / r0n;
    REAL alpha = -2.0 * energy / mu;

    if (dt == (REAL)0.0) return;

    REAL chi;
    if (std::abs(alpha) < 1e-12) {
        chi = std::sqrt(mu) * dt / r0n;
    } else if (alpha > 0.0) {
        chi = std::sqrt(mu) * dt * alpha;
    } else {
        REAL sign_dt = (dt >= 0.0 ? 1.0 : -1.0);
        REAL a = 1.0 / alpha; 
        REAL sqrtm = std::sqrt(-a);
        REAL arg = -2.0 * mu * alpha * dt / (vr0 + sign_dt * std::sqrt(-mu / alpha) * (1.0 - r0n * alpha));
        if (arg > 0.0)
            chi = sign_dt * sqrtm * std::log(arg);
        else
            chi = sign_dt * std::sqrt(mu) * std::abs(alpha) * dt; // fallback
    }

    auto eval_F_and_Fprime = [&](REAL chi_in, REAL &F_out, REAL &Fp_out) {
        REAL z = alpha * chi_in * chi_in;
        REAL C = stumpff_C(z);
        REAL S = stumpff_S(z);
        REAL sqrtmu = std::sqrt(mu);
        REAL p = r0n * vr0 / sqrtmu;
        REAL q = 1.0 - alpha * r0n;

        F_out = p * chi_in*chi_in * C + q * chi_in*chi_in*chi_in * S + r0n * chi_in - sqrtmu * dt;
        Fp_out = p * chi_in * (1.0 - alpha * chi_in*chi_in * S) + q * chi_in*chi_in * C + r0n;
    };

    REAL correction = 1.0;
    int iter = 0;
    for (; iter < maxiter; ++iter) {
        REAL F, Fp;
        eval_F_and_Fprime(chi, F, Fp);

        REAL h = std::max((REAL)1e-6, (REAL)1e-8 * std::max((REAL)1.0, std::abs(chi)));

        REAL Fp_ph, Fp_mh, F_ph, F_mh;
        eval_F_and_Fprime(chi + h, F_ph, Fp_ph);
        eval_F_and_Fprime(chi - h, F_mh, Fp_mh);

        REAL Fpp = (Fp_ph - Fp_mh) / (2.0 * h);
        REAL Fppp = (Fp_ph - 2.0*Fp + Fp_mh) / (h*h);

        if (std::abs(Fp) < 1e-20) { 
            correction = -F * 1e20;
            chi += correction;
            break;
        }

        REAL delta1 = - F / Fp;
        REAL denom2 = Fp + 0.5 * delta1 * Fpp;
        REAL delta2 = (std::abs(denom2) < 1e-20) ? delta1 : - F / denom2;
        REAL denom3 = Fp + 0.5 * delta1 * (Fpp + (delta1*delta1 * Fppp) / 3.0);
        REAL delta3 = (std::abs(denom3) < 1e-20) ? delta2 : - F / denom3;

        correction = delta3;
        chi += correction;

        if (std::abs(correction) < tol) break;
    }

    REAL z = alpha * chi * chi;
    REAL C = stumpff_C(z);
    REAL S = stumpff_S(z);
    REAL sqrtmu = std::sqrt(mu);

    REAL f_gauss = 1.0 - (chi*chi / r0n) * C;
    REAL g_gauss = dt - (chi*chi*chi / sqrtmu) * S;

    vec3 r = f_gauss * r0 + g_gauss * v0;
    REAL rn = r.magnitude();

    REAL fdot = (sqrtmu / (rn * r0n)) * (alpha * chi*chi*chi * S - chi);
    REAL gdot = 1.0 - (chi*chi / rn) * C;

    vec3 v = fdot * r0 + gdot * v0;

    r0 = r;
    v0 = v;
}

void update_WH_planetary(Simulation* _s) {
    const size_t N = _s->bodies.positions.size();
    if (N <= 1) return;

    const REAL G = _s->get_G();
    const REAL dt = _s->get_dt();
    _s -> solver -> set_bounds(1, -1);
    vec3 r_central = _s->bodies.positions[0];
    vec3 v_central = _s->bodies.velocities[0];
    REAL m_central = _s->bodies.masses[0];

    parallel_for(size_t(1), N, [_s, r_central, v_central, m_central, G, dt](size_t i){
        vec3 r_rel = _s->bodies.positions[i] - r_central;
        vec3 v_rel = _s->bodies.velocities[i] - v_central;

        REAL mu_i = G * (m_central + _s->bodies.masses[i]);
        update_kepler(r_rel, v_rel, mu_i, dt*0.5);
        _s->bodies.positions[i] = r_central + r_rel;
        _s->bodies.velocities[i] = v_central + v_rel;
    });

    _s -> solver -> compute_forces();

    parallel_for(size_t(1), N, [_s, r_central, v_central, m_central, G, dt](size_t i){
        _s->bodies.velocities[i] += _s->bodies.accs[i] * dt;
        vec3 r_rel = _s->bodies.positions[i] - r_central;
        vec3 v_rel = _s->bodies.velocities[i] - v_central;

        REAL mu_i = G * (m_central + _s->bodies.masses[i]);
        update_kepler(r_rel, v_rel, mu_i, dt*0.5);

        _s->bodies.positions[i] = r_central + r_rel;
        _s->bodies.velocities[i] = v_central + v_rel;
        _s->bodies.accs[i].reset(); 
    });
}

void Solver::set_force(force_type _t){
    auto idx = static_cast<int>(_t);
    if (critical_if(idx > 3 || idx < 0, "Invalid force type")) return;
    

    switch (_t){
        case NEWTON:
            this -> get_force = newtonian;
        case PN1:
            this -> get_force = pn1;
        case PN2:
            this -> get_force = pn2;
        case PN25:
            this -> get_force = pn25;
    }
}


inline void get_new_temp(Simulation* _s, size_t body, vec3 pos, vec3 vel, REAL temp, REAL mass){ 
    /*
    * we are trying to compute the delta t so we use the formula 
    *  21 * G * M * m * n * e²
    *  ------------------------ = E_t
    *         2 * r^6

    * from here it is trival, we can just divide by the mass and a coeff and multiply 
    * byt dt to get the delta t 
    */
    // sketchy edge case...
    if (!_s -> bodies.get_acc_of(body).sqr_magn()) return;

    vec3 rel_v = _s -> bodies.get_velocity_of(body) - vel;
    REAL relv_sq = rel_v.sqr_magn();
    REAL r = (_s -> bodies.get_position_of(body) - pos).magnitude();
    
    /*
    * omega is the angular velocity of the body
    * thus it can be calculated by getting the perpendicular
    * vector to the projection of the velocity onto the acceleration
    * assuming the acceleration comes purely from the mutual
    * attraction between the bodies
    * ergo 
    * W  = V x a / a x a * vec(a)
    */
    REAL w_squared = (_s -> bodies.get_acc_of(body) * (rel_v * _s -> bodies.get_acc_of(body)) / (_s -> bodies.get_acc_of(body).sqr_magn())).sqr_magn();

    REAL eccentricity = get_eccentricity(_s, body, relv_sq, w_squared, r, mass);

    REAL delta_temperature = 21.0 / 2.0 * _s -> get_G() * mass * _s -> bodies.masses[body] * eccentricity * eccentricity * w_squared / (r*r*r*r*r*r + _s->softening);
    delta_temperature *= _s -> get_dt(); // we want it to be per unit of time
    delta_temperature /= _s -> bodies.get_mass_of(body) * _s -> get_heat_capacity(); // then we get how many K the body got in dt 

    // boltzmann
    delta_temperature -= _s -> get_boltzmann() * std::pow(_s -> bodies.get_temp_of(body), 4) / ( _s -> get_heat_capacity() * _s -> bodies.get_mass_of(body));

    _s -> bodies.get_temp_of(body) += delta_temperature; // done! 
}


/*
 
void adaptive_euler(Simulation* _s){
    for (auto& body : _s -> bodies){
        REAL magn = body.acceleration.magnitude();
        T r = body.acceleration / magn;

        REAL E_i = _s -> bodies.get_velocity_of(body).sqr_magn() * .5 *body.mass - std::sqrt(_s -> get_G() * (_s -> get_total_mass() - body.mass) / magn);

        Body rollback = body;
        update_euler(_s);

        REAL E_f = _s -> bodies.get_velocity_of(body).sqr_magn() * .5 *body.mass - std::sqrt(_s -> get_G() * (_s -> get_total_mass() - body.mass) / body.acceleration.magnitude());
        REAL iters = std::log10(std::abs((E_i - E_f) / E_i));

        if (iters < _s -> get_adaptive_coeff())
            return;

        iters = _s -> get_adaptive_coeff() - iters;

        _s -> bodies.get_position_of(body) = rollback.position;
        body.acceleration = rollback.acceleration;
        body.temp = rollback.temp;
        _s -> bodies.get_velocity_of(body) = rollback.velocity;

        _s -> set_dt(_s -> get_dt() / iters);
        for (REAL i = 0; i < iters; ++i){
            update_euler(_s);
        }

        _s -> set_dt(_s -> get_dt() * iters);
    }
}*/
//===---------------------------------------------------------===//
// Force calculation methods                                     //
//===---------------------------------------------------------===//


vec3 newtonian(REAL m1, REAL m2, vec3 v1, vec3 v2, vec3 p1, vec3 p2, Simulation* _s){
    vec3 n = (p2 - p1);
    REAL r = n.magnitude();

    return n.normalize() * _s -> get_G() *m1*m2/(r*r+_s -> softening);
}

//===---------------------------------------------------------===//
// Force calculation methods                                     //
//===---------------------------------------------------------===//

inline std::vector<vec3> get_new_pos(vec3* position, vec3* velocity, vec3* acceleration, REAL step){
    std::vector<vec3> retval;
    retval.push_back(velocity -> update_by(acceleration, step));
    retval.push_back(position -> update_by(&retval[0], step));
    return retval;
}



}