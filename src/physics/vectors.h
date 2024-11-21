#pragma once

struct vec2 {
    double x;
    double y;

    vec2() : x(0), y(0) {}

    vec2(double x, double y) : x(x), y(y) {}

    vec2 operator*(double scalar) const;
    vec2 operator/(double scalar) const;
    vec2 operator/(float scalar) const;
    vec2 operator/(const int& other) const;
    vec2 operator+(const vec2& other) const;
    vec2 operator-(const vec2& other) const;
    vec2 operator-() const;
    vec2& operator+=(const vec2& other);
    bool operator==(const vec2& other) const;
    vec2 operator*(const vec2& other) const;
    vec2 operator*(const int& other) const;
    float magnitude();
    vec2 direction(vec2 v2);
    vec2 normalize();

};

#include "vectors_helper.h"
