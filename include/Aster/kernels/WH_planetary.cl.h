#pragma once
#include <string>

namespace Aster{
namespace GPU{

inline std::string wh_cl_3d = R"CLC(
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#define tol 1e-16
#define maxiter 500
#define z_series_thresh 1e-8

static double stumpff_C(double z){
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

static double stumpff_S(double z){
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
static void eval_F_and_Fprime(double chi, double alpha, double p, double q, double r0n, double sqrtmu, double dt, double *F_out, double *Fp_out) {
    double z = alpha * chi * chi;
    double C = stumpff_C(z);
    double S = stumpff_S(z);

    double F = p * chi*chi * C + q * chi*chi*chi * S + r0n * chi - sqrtmu * dt;
    double Fp = p * chi * (1.0 - alpha * chi*chi * S) + q * chi*chi * C + r0n;

    *F_out = F;
    *Fp_out = Fp;
}

static void update_kepler(double3 *r0_ptr, double3 *v0_ptr, double mu, double dt) {
    double3 r0 = *r0_ptr;
    double3 v0 = *v0_ptr;

    double r0n = length(r0);
    if (r0n == 0.0) {
        return;
    }
    double v0n = length(v0);
    double vr0 = dot(r0, v0) / r0n;

    double energy = v0n*v0n/2.0 - mu / r0n;
    double alpha = -2.0 * energy / mu;

    if (dt == 0.0) {
        return;
    }

    double chi = 0.0;
    if (fabs(alpha) < 1e-12) {
        chi = sqrt(mu) * dt / r0n;
    } else if (alpha > 0.0) {
        chi = sqrt(mu) * dt * alpha;
    } else {
        double sign_dt = (dt >= 0.0 ? 1.0 : -1.0);
        double a = 1.0 / alpha;
        double sqrtm = sqrt(-a);
        double denom = (vr0 + sign_dt * sqrt(-mu / alpha) * (1.0 - r0n * alpha));
        double arg = 0.0;
        if (denom != 0.0) arg = -2.0 * mu * alpha * dt / denom;
        if (arg > 0.0)
            chi = sign_dt * sqrtm * log(arg);
        else
            chi = sign_dt * sqrt(mu) * fabs(alpha) * dt; 
    }

    double sqrtmu = sqrt(mu);
    double p = r0n * vr0 / sqrtmu;
    double q = 1.0 - alpha * r0n;

    double correction = 1.0;
    int iter = 0;
    for (; iter < maxiter; ++iter) {
        double F, Fp;
        eval_F_and_Fprime(chi, alpha, p, q, r0n, sqrtmu, dt, &F, &Fp);

        double h = fmax(1e-6, 1e-8 * fmax(1.0, fabs(chi)));
        double F_ph, Fm_ph;
        double Fp_ph, Fp_mh;
        eval_F_and_Fprime(chi + h, alpha, p, q, r0n, sqrtmu, dt, &F_ph, &Fp_ph);
        eval_F_and_Fprime(chi - h, alpha, p, q, r0n, sqrtmu, dt, &Fm_ph, &Fp_mh);

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

    double f_gauss = 1.0 - (chi*chi / r0n) * C;
    double g_gauss = dt - (chi*chi*chi / sqrtmu) * S;

    double3 r = f_gauss * r0 + g_gauss * v0;
    double rn = length(r);

    double fdot = (sqrtmu / (rn * r0n)) * (alpha * chi*chi*chi * S - chi);
    double gdot = 1.0 - (chi*chi / rn) * C;

    double3 v = fdot * r0 + gdot * v0;

    *r0_ptr = r;
    *v0_ptr = v;
}

__kernel void wh_planetary(
    __global double* pos, 
    __global double* vel,
    __global double* acc,
    __global double* masses,
    const double dt,
    const double G,
    const int lower,
    const int upper,
    const int pass
) {
    const uint gid = get_global_id(0);
    if ((int)gid < lower) return;
    if ((int)gid >= upper) return;

    double3 r_central = vload3(0, pos);
    double3 v_central = vload3(0, vel);
    double m_central = masses[0];

    double3 p = vload3(gid, pos);
    double3 v = vload3(gid, vel);
    double3 a = vload3(gid, acc);

    if (pass) {
        v += a * dt;
    }

    double3 r_rel = p - r_central;
    double3 v_rel = v - v_central;

    double mu_i = G * (m_central + masses[gid]);

    update_kepler(&r_rel, &v_rel, mu_i, dt * 0.5);

    p = r_central + r_rel;
    v = v_central + v_rel;

    vstore3(p, gid, pos);
    vstore3(v, gid, vel);
    vstore3((double3)(0.0,0.0,0.0), gid, acc);
}
    
)CLC";



}
}
