#pragma once

namespace Aster{
namespace Units{
constexpr double c    = 299792458;
constexpr double G    = 6.6743015e-11;

  //============================================//
 //                UNITS OF TIME               //
//============================================//

constexpr double second = 1.0;
constexpr double minute = 60.0 * second;
constexpr double hour   = 60.0 * minute;
constexpr double day    = 24.0 * hour;
constexpr double year   = 365.25 * day;
constexpr double kyr    = 1e3 * year;
constexpr double Myr    = 1e6 * year;
constexpr double Gyr    = 1e9 * year;

constexpr double operator"" _s(long double v){ return static_cast<double>(v) * second; }
constexpr double operator"" _min(long double v){ return static_cast<double>(v) * minute; }
constexpr double operator"" _h(long double v){ return static_cast<double>(v) * hour; }
constexpr double operator"" _d(long double v){ return static_cast<double>(v) * day; }
constexpr double operator"" _yr(long double v){ return static_cast<double>(v) * year; }
constexpr double operator"" _kyr(long double v){ return static_cast<double>(v) * kyr; }
constexpr double operator"" _Myr(long double v){ return static_cast<double>(v) * Myr; }
constexpr double operator"" _Gyr(long double v){ return static_cast<double>(v) * Gyr; }
constexpr double operator"" _yr(unsigned long long v){ return static_cast<double>(v) * year; }
constexpr double operator"" _Myr(unsigned long long v){ return static_cast<double>(v) * Myr; }
constexpr double operator"" _Gyr(unsigned long long v){ return static_cast<double>(v) * Gyr; }

  //============================================//
 //                UNITS OF MASS               //
//============================================//

constexpr double kilogram = 1.0;
constexpr double gram     = 1e-3 * kilogram;
constexpr double solar_mass = 1.98847e30 * kilogram;
constexpr double earth_mass = 5.9722e24 * kilogram;
constexpr double jupiter_mass = 1.89813e27 * kilogram;

constexpr double operator"" _kg(long double v){ return static_cast<double>(v) * kilogram; }
constexpr double operator"" _g(long double v){ return static_cast<double>(v) * gram; }
constexpr double operator"" _Msun(long double v){ return static_cast<double>(v) * solar_mass; }
constexpr double operator"" _Mearth(long double v){ return static_cast<double>(v) * earth_mass; }
constexpr double operator"" _Mjup(long double v){ return static_cast<double>(v) * jupiter_mass; }
constexpr double operator"" _Msun(unsigned long long v){ return static_cast<double>(v) * solar_mass; }

  //============================================//
 //                UNITS OF SPACE              //
//============================================//
 
constexpr double meter = 1.0;

constexpr double millimeter = 1e-3 * meter;
constexpr double centimeter = 1e-2 * meter;
constexpr double kilometer  = 1e3  * meter;
constexpr double astronomical_unit = 1.495978707e11 * meter; // AU (IAU 2012)
constexpr double light_second      = c * second;
constexpr double light_minute      = c * minute;
constexpr double light_hour        = c * hour;
constexpr double light_day         = c * day;
constexpr double light_year        = c * year;

constexpr double parsec            = 3.0856775814913673e16 * meter;
constexpr double kiloparsec        = 1e3 * parsec;
constexpr double megaparsec        = 1e6 * parsec;

constexpr double operator"" _m(long double v) {    return static_cast<double>(v) * meter;}
constexpr double operator"" _mm(long double v) {    return static_cast<double>(v) * millimeter;}
constexpr double operator"" _cm(long double v) {    return static_cast<double>(v) * centimeter;}
constexpr double operator"" _km(long double v) {    return static_cast<double>(v) * kilometer;}
constexpr double operator"" _AU(long double v) {    return static_cast<double>(v) * astronomical_unit;}
constexpr double operator"" _ls(long double v) {    return static_cast<double>(v) * light_second;}
constexpr double operator"" _ld(long double v) {    return static_cast<double>(v) * light_day;}
constexpr double operator"" _ly(long double v) {    return static_cast<double>(v) * light_year;}
constexpr double operator"" _pc(long double v) {    return static_cast<double>(v) * parsec;}
constexpr double operator"" _kpc(long double v) {    return static_cast<double>(v) * kiloparsec;}
constexpr double operator"" _Mpc(long double v) {    return static_cast<double>(v) * megaparsec;}



constexpr double H0 = 70.0 * kilometer / second / megaparsec;

constexpr double pi = 3.14159265358979323846;
constexpr double two_pi = 2.0 * pi;
constexpr double half_pi = 0.5 * pi;

constexpr double hubble_time = 1.0 / H0;
constexpr double critical_density = 3.0 * H0 * H0 / (8.0 * pi * G);

constexpr double machine_epsilon = 2.2204460492503131e-16;

}
}