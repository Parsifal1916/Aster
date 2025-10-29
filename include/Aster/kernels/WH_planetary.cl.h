#pragma once
#include <string>

namespace Aster{
namespace GPU{

inline std::string wh_cl_3d = R"CLC(
#pragma OPENCL EXTENSION cl_khr_fp64 : enable


inline double stumpff_C(double z) {
    const double z_series_thresh = 1e-8;
    if (fabs(z) < z_series_thresh) {
        double z2 = z*z;
        double z3 = z2*z;
        return 0.5 - z/24.0 + z2/720.0 - z3/40320.0;
    } else if (z > 0.0) {
        double s = sqrt(z);
        return (1.0 - cos(s)) / z;
    } else {
        double s = sqrt(-z);
        return (cosh(s) - 1.0) / (-z);
    }
}

inline double stumpff_S(double z) {
    const double z_series_thresh = 1e-8;
    if (fabs(z) < z_series_thresh) {
        double z2 = z*z;
        double z3 = z2*z;
        return (1.0/6.0) - z/120.0 + z2/5040.0 - z3/362880.0;
    } else if (z > 0.0) {
        double s = sqrt(z);
        return (s - sin(s)) / (s*s*s);
    } else {
        double s = sqrt(-z);
        return (sinh(s) - s) / (s*s*s);
    }
}


inline void eval_F(double chi_in, double *F, double *Fp,
                   double r0n, double vr0, double alpha,
                   double mu, double dt) {
    double z = alpha * chi_in * chi_in;
    double C = stumpff_C(z);
    double S = stumpff_S(z);
    double sqrtmu = sqrt(mu);
    double p = r0n * vr0 / sqrtmu;
    double q = 1.0 - alpha * r0n;

    *F  = p * chi_in*chi_in * C + q * chi_in*chi_in*chi_in * S + r0n * chi_in - sqrtmu * dt;
    *Fp = p * chi_in * (1.0 - alpha * chi_in*chi_in * S) + q * chi_in*chi_in * C + r0n;
}


inline double3 kepler_update(double3 r0, double3 v0, double mu, double dt, __private double3* v_out) {
    const double tol_base = 1e-16;
    const int maxiter = 500;

    double r0n = length(r0);
    double v0n = length(v0);

    double vr0 = dot(r0, v0) / r0n;

    double energy = 0.5 * v0n * v0n - mu / r0n;
    double alpha = -2.0 * energy / mu;

    if (dt == 0.0) {
        *v_out = v0;
        return r0;
    }

    double chi;
    if (fabs(alpha) < 1e-12) {
        chi = sqrt(mu) * dt / r0n;
    } else if (alpha > 0.0) {
        chi = sqrt(mu) * dt * alpha;
    } else {
        double sign_dt = (dt >= 0.0 ? 1.0 : -1.0);
        double a = 1.0 / alpha; 
        double sqrtm = sqrt(-a);
        double arg = -2.0 * mu * alpha * dt / (vr0 + sign_dt * sqrt(-mu / alpha) * (1.0 - r0n * alpha));
        if (arg > 0.0)
            chi = sign_dt * sqrtm * log(arg);
        else
            chi = sign_dt * sqrt(mu) * fabs(alpha) * dt;
    }

    double tol = tol_base;

    double correction = 1.0;
    for (int iter = 0; iter < maxiter; ++iter) {
        double F, Fp;
        eval_F(chi, &F, &Fp, r0n, vr0, alpha, mu, dt);

       double h = fmax((double)1e-6, (double)1e-8 * fmax((double)1.0, fabs(chi)));

        double F_ph, F_mh, Fp_ph, Fp_mh;
        eval_F(chi + h, &F_ph, &Fp_ph, r0n, vr0, alpha, mu, dt);
        eval_F(chi - h, &F_mh, &Fp_mh, r0n, vr0, alpha, mu, dt);

        double Fpp = (Fp_ph - Fp_mh) / (2.0 * h);
        double Fppp = (Fp_ph - 2.0*Fp + Fp_mh) / (h*h);

        if (fabs(Fp) < 1e-20) {
            correction = -F * 1e20;
            chi += correction;
            break;
        }

        double delta1 = - F / Fp;
        double denom2 = Fp + 0.5 * delta1 * Fpp;
        double delta2 = (fabs(denom2) < 1e-20) ? delta1 : - F / denom2;
        double denom3 = Fp + 0.5 * delta1 * (Fpp + (delta1*delta1 * Fppp) / 3.0);
        double delta3 = (fabs(denom3) < 1e-20) ? delta2 : - F / denom3;

        correction = delta3;
        chi += correction;

        if (fabs(correction) < tol) break;
    }

    double z = alpha * chi * chi;
    double C = stumpff_C(z);
    double S = stumpff_S(z);
    double sqrtmu = sqrt(mu);

    double chi2 = chi * chi;
    double chi3 = chi2 * chi;

    double f_gauss = 1.0 - (chi2 / r0n) * C;
    double g_gauss = dt - (chi3 / sqrtmu) * S;

    double3 r = f_gauss * r0 + g_gauss * v0;
    double rn = length(r);

    double fdot = (sqrtmu / (rn * r0n)) * (alpha * chi3 * S - chi);
    double gdot = 1.0 - (chi2 / rn) * C;

    double3 v = fdot * r0 + gdot * v0;

    *v_out = v;
    return r;
}

__kernel void wh_first(
    const uint N,
    const double G,
    const double dt,
    __global const double* masses,
    __global double* positions,
    __global double* velocities
){
    uint i = get_global_id(0);
    i++;
    if (i >= N) return;

    double3 r_central = vload3(0, positions);
    double3 v_central = vload3(0, velocities); 
    double m_central = masses[0];   

    double3 r = vload3(i, positions) - r_central;
    double3 v = vload3(i, velocities) - v_central;

    double mu_i = G * (m_central + masses[i]);
    double3 v_new;
    double3 r_new = kepler_update(r, v, mu_i, dt*0.5, &v_new);

    vstore3(r_central + r_new, i, positions);
    vstore3(v_central + v_new, i, velocities);
}

__kernel void wh_second(
    const uint N,
    const double G,
    const double dt,
    __global const double* masses,
    __global double* positions,
    __global double* velocities,
    __global double* accs
){
    uint i = get_global_id(0);
    i++;
    if (i >= N) return;

    double3 r_central = vload3(0, positions);
    double3 v_central = vload3(0, velocities); 
    double m_central = masses[0];  

    double3 r = vload3(i, positions) - r_central;
    double3 v = vload3(i, velocities) - v_central;    
    double3 a = vload3(i, accs);
    v += a * dt; 

    double mu_i = G * (m_central + masses[i]);
    double3 v_new;
    double3 r_new = kepler_update(r, v, mu_i, dt*0.5, &v_new);

    vstore3(r_central + r_new, i, positions);
    vstore3(v_central + v_new, i, velocities);
    vstore3((double3)(0.0,0.0,0.0), i, accs);
}
    
)CLC";



}
}
