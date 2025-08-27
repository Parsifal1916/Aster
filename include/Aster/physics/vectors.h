#pragma once

#include <iostream>

using REAL = double;

namespace Aster{

struct vec2 {
    REAL x;
    REAL y;

    vec2() : x((REAL)0), y((REAL)0) {}

    vec2(REAL scalar) : x(scalar), y(scalar) {}
    vec2(REAL x, REAL y) : x(x), y(y) {}

    vec2 operator*(double scalar) const;
    vec2 operator/(double scalar) const;
    vec2 operator/(float scalar) const;
    vec2 operator/(const int& other) const;
    vec2 operator+(const vec2& other) const;
    vec2 operator-(const vec2& other) const;
    vec2 operator-() const;
    REAL& operator[](size_t index);
    vec2& operator+=(const vec2& other);
    bool operator==(const vec2& other) const;
    REAL operator*(const vec2& other) const;
    vec2 operator*(const int& other) const;
    vec2& operator/=(const REAL& other);
    REAL magnitude() const;
    vec2 direction(vec2 v2);
    vec2 normalize();
    bool is_fine() const;
    REAL sqr_magn() const;
    void reset();
    vec2 update_by(vec2* v, REAL delta);

};



struct vec3 {
    REAL x;
    REAL y;
    REAL z;

    vec3() : x(0), y(0), z(0) {}

    vec3(REAL scalar) : x(scalar), y(scalar), z(scalar) {}
    vec3(REAL x, REAL y, REAL z) : x(x), y(y), z(z) {}

    vec3 operator*(double scalar) const;
    vec3 operator/(double scalar) const;
    vec3 operator/(float scalar) const;
    vec3 operator/(const int& other) const;
    vec3 operator+(const vec3& other) const;
    vec3 operator-(const vec3& other) const;
    vec3 operator-() const;
    REAL& operator[](size_t index);
    vec3& operator+=(const vec3& other);
    bool operator==(const vec3& other) const;
    REAL operator*(const vec3& other) const;
    vec3 operator*(const int& other) const;
    vec3& operator/=(const REAL& other);
    REAL magnitude() const;
    bool is_fine() const;
    vec3 update_by(vec3* v, REAL delta);
    vec3 direction(vec3 v3);
    vec3 normalize();
    REAL sqr_magn() const;
    void reset();
};

vec2 operator*(REAL scalar, vec2 v);

vec3 operator*(REAL scalar, vec3 v);

std::ostream& operator<<(std::ostream& os, const vec3& v);
std::ostream& operator<<(std::ostream& os, const vec2& v);
}