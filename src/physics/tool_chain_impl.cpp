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

 
inline vec3 get_pn25(REAL m1, REAL m2, vec3 v1, vec3 v2, vec3 p1, vec3 p2, Simulation* _s, REAL inv_r, vec3 norm){
    const REAL G = _s -> get_G();
    REAL pseudo_newton1 = G * m1 * inv_r;
    REAL pseudo_newton2 = G * m2 * inv_r;

    REAL sq_v_diff = (v1 - v2).sqr_magn();

    vec3 internal_1 = (v1-v2) * (-sq_v_diff + 2 * pseudo_newton1 - 8 * pseudo_newton2);

    vec3 internal_2 = norm * (norm*v1 - norm*v2) * (3* sq_v_diff - 6 * pseudo_newton1 + 52.0/3.0 * pseudo_newton2);

    return (internal_1 + internal_2) * 4.0/5.0 * G *G *m1*m2 * inv_r*inv_r*inv_r ;
}

 
inline vec3 get_pn2(REAL m1, REAL m2, vec3 v1, vec3 v2, vec3 p1, vec3 p2, Simulation* _s, REAL inv_r, vec3 norm){
    const REAL v2_sq = v2.sqr_magn();
    const REAL v1_sq = v1.sqr_magn();

    const REAL v_prod = v1*v2;
    const REAL G = _s -> get_G();

    REAL nv1 = norm * v1;
    REAL nv2 = norm * v2;
    REAL nv22 = nv2 * nv2;
    REAL pseudo_newtonian1 = G * m1 * inv_r;
    REAL pseudo_newtonian2 = G * m2 * inv_r;
    const REAL nv12 = nv1 * nv1;

    return G * m2 * inv_r * inv_r * (norm * (-2 *v2_sq * v2_sq + 4 * v2_sq * v_prod - 2 *v_prod * v_prod
    + 3.0/2.0 * v1_sq * nv22 + 9.0/2.0 * v2_sq * nv22 - 6 * v_prod * nv22
    -15.0/8.0 * nv22 * nv22 + pseudo_newtonian1 * (-15.0/4.0 * v1_sq+5.0/4.0 *v2_sq-5.0/2.0 * v_prod
    +39.0/2.0 * nv12-39.0 * nv1 * nv2+ 17.0/2.0 * nv22)
    + pseudo_newtonian2 * (4 * v2_sq- 8 * v_prod +2* nv12
    -4 * nv1 * nv2 - 6 * nv22))
    +(v1-v2) * (v1_sq * nv2+ 4 * v2_sq * nv1- 5 * v2_sq * nv2
    -4 * v_prod * nv1+4 * v_prod * nv2 - 6 * nv1 * nv22
    +9.0/2.0 * nv22 * nv2 + pseudo_newtonian1 * (-63.0/4.0 * nv1+55.0/4.0 * nv2)
    + pseudo_newtonian2 * (-2 * nv1- 2 * nv2)))
    +G*G*G * m2 * inv_r*inv_r*inv_r*inv_r * norm * (-57.0/4.0 *  m1*m1 - 9 * m2 * m2 - 69.0/2.0 * m1 * m2);
}

 
inline vec3 get_pn1(REAL m1, REAL m2, vec3 v1, vec3 v2, vec3 p1, vec3 p2, Simulation* _s, REAL inv_r, vec3 norm){
    const REAL v2_sq = v2.sqr_magn();
    const REAL G = _s -> get_G(); 
    vec3 a = norm * (-v1.sqr_magn() - 2 * v2_sq + 4 *v1*v2 + 3.0/2.0 * ((norm *v2) * (norm * v2)) + 5 * G * m1 * inv_r + 4 * G * m2 *inv_r);
    a += (v1-v2) * (4.0 * norm * v1 - 3 * norm * v2);
    a = a * G * m2 * inv_r * inv_r;

    return a;
}

vec3 compute_gravitational(REAL m1, REAL m2, vec3 v1, vec3 v2, vec3 p1, vec3 p2, Simulation* _s){
    vec3 d = p2 - p1;
    REAL r2 = d.x*d.x + d.y*d.y + d.z*d.z + _s->softening;
    REAL inv_r = 1.0/std::sqrt(r2);
    REAL inv_r2 = inv_r * inv_r;
    vec3 norm = d * inv_r;

    vec3 newton = norm * (_s -> get_G() * m2 * m1 * inv_r2);
    REAL inv_c = 1.0 / _s -> get_c();

    vec3 a1(0), a2(0), a25(0);
    auto index = _s->solver->force_law;

    if (index.pn1())  a1  = inv_c * inv_c * get_pn1(m1, m2, v1, v2, p1, p2, _s, inv_r, -norm);
    if (index.pn2())  a2  = inv_c * inv_c * inv_c * inv_c * get_pn2(m1, m2, v1, v2, p1, p2, _s, inv_r, -norm);
    if (index.pn25()) a25 = inv_c * inv_c * inv_c * inv_c * inv_c * get_pn25(m1, m2, v1, v2, p1, p2, _s, inv_r, -norm);
    return newton + a1 * m1 + a2 * m1 + a25 * m1;
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
    if (dt == 0.0) return;

    const REAL tol = 1e-15;
    const int maxiter = 100;
    const REAL r0n = r0.magnitude();
    const REAL v0n = v0.magnitude();
    const REAL rdotv = r0.x * v0.x + r0.y * v0.y + r0.z * v0.z;
    const REAL sqrtmu = std::sqrt(mu);
    
    const REAL energy = 0.5 * v0n * v0n - mu / r0n;
    const REAL alpha = -2.0 * energy / mu;

    auto get_stumpff = [](REAL z, REAL& c, REAL& s) {
        if (z > 1e-6) {
            REAL sqrtz = std::sqrt(z);
            c = (1.0 - std::cos(sqrtz)) / z;
            s = (sqrtz - std::sin(sqrtz)) / (sqrtz * z);
        } else if (z < -1e-6) {
            REAL sqrtz = std::sqrt(-z);
            c = (std::cosh(sqrtz) - 1.0) / (-z);
            s = (std::sinh(sqrtz) - sqrtz) / ((-z) * sqrtz);
        } else {
            c = 1.0 / 2.0 - z / 24.0 + z * z / 720.0;
            s = 1.0 / 6.0 - z / 120.0 + z * z / 5040.0;
        }
    };

    REAL chi = (std::abs(alpha) < 1e-10) ? (sqrtmu * dt / r0n) : (sqrtmu * dt * alpha);
    
    for (int i = 0; i < maxiter; ++i) {
        REAL z = alpha * chi * chi;
        REAL c, s;
        get_stumpff(z, c, s);

        REAL f = r0n * chi * (1.0 - z * s) + rdotv / sqrtmu * chi * chi * c + chi - sqrtmu * dt;
        REAL fp = r0n * (1.0 - z * c) + rdotv / sqrtmu * chi * (1.0 - z * s) + 1.0; 
        
        REAL f_val = r0n * chi * (1.0 - alpha * chi * chi * s) + (rdotv / sqrtmu) * chi * chi * c + r0n * chi - sqrtmu * dt; // Wait, standard form:
        
        REAL q = 1.0 - alpha * r0n;
        REAL p = rdotv / sqrtmu;
        
        REAL current_F = p * chi * chi * c + q * chi * chi * chi * s + r0n * chi - sqrtmu * dt;
        REAL current_Fp = p * chi * (1.0 - z * s) + q * chi * chi * c + r0n;
        REAL current_Fpp = p * (1.0 - z * c) + q * chi * (1.0 - z * s);

        REAL delta = current_F / (current_Fp - 0.5 * current_F * current_Fpp / current_Fp);
        chi -= delta;
        if (std::abs(delta) < tol) break;
    }

    REAL z = alpha * chi * chi;
    REAL c, s;
    get_stumpff(z, c, s);

    REAL f = 1.0 - (chi * chi / r0n) * c;
    REAL g = dt - (chi * chi * chi / sqrtmu) * s;

    vec3 r = f * r0 + g * v0;
    REAL rn = r.magnitude();

    REAL fdot = (sqrtmu / (rn * r0n)) * (alpha * chi * chi * chi * s - chi);
    REAL gdot = 1.0 - (chi * chi / rn) * c;

    v0 = fdot * r0 + gdot * v0;
    r0 = r;
}

void update_WH_planetary(Simulation* _s) {
    const size_t N = _s->bodies.positions.size();
    if (N <= 1) return;

    const REAL G = _s->get_G();
    const REAL dt = _s->get_dt();
    const REAL mu = G * _s->bodies.masses[0];

    parallel_for(size_t(1), N, [_s, mu, dt](size_t i) {
        vec3 r_rel = _s->bodies.positions[i] - _s->bodies.positions[0];
        vec3 v_rel = _s->bodies.velocities[i] - _s->bodies.velocities[0];
        update_kepler(r_rel, v_rel, mu, dt * 0.5);
        _s->bodies.positions[i] = _s->bodies.positions[0] + r_rel;
        _s->bodies.velocities[i] = _s->bodies.velocities[0] + v_rel;
    });

    _s->bodies.positions[0] += _s->bodies.velocities[0] * (dt * 0.5);

    _s->solver->compute_forces();

    parallel_for(size_t(0), N, [_s, dt](size_t i) {
        _s->bodies.velocities[i] += _s->bodies.accs[i] * dt;
    });

    _s->bodies.positions[0] += _s->bodies.velocities[0] * (dt * 0.5);

    parallel_for(size_t(1), N, [_s, mu, dt](size_t i) {
        vec3 r_rel = _s->bodies.positions[i] - _s->bodies.positions[0];
        vec3 v_rel = _s->bodies.velocities[i] - _s->bodies.velocities[0];
        update_kepler(r_rel, v_rel, mu, dt * 0.5);
        _s->bodies.positions[i] = _s->bodies.positions[0] + r_rel;
        _s->bodies.velocities[i] = _s->bodies.velocities[0] + v_rel;
    });

    parallel_for(size_t(0), N, [_s](size_t i) {
        _s->bodies.accs[i].reset();
    });
}

void Solver::set_force(PN_expansion_type _t){
    this->force_law = _t;
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