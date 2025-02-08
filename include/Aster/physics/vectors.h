#pragma once

namespace Aster{

struct vec2 {
    double x;
    double y;

    vec2() : x(0), y(0) {}

    vec2(double scalar) : x(scalar), y(scalar) {}
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
    double operator*(const vec2& other) const;
    vec2 operator*(const int& other) const;
    vec2& operator/=(const double& other);
    float magnitude();
    vec2 direction(vec2 v2);
    vec2 normalize();
    double sqr_magn() const;
    void reset();
    vec2 update_by(vec2* v, double delta);

};

struct vec3 {
    double x;
    double y;
    double z;

    vec3() : x(0), y(0), z(0) {}

    vec3(double scalar) : x(scalar), y(scalar), z(scalar) {}
    vec3(double x, double y, double z) : x(x), y(y), z(z) {}

    vec3 operator*(double scalar) const;
    vec3 operator/(double scalar) const;
    vec3 operator/(float scalar) const;
    vec3 operator/(const int& other) const;
    vec3 operator+(const vec3& other) const;
    vec3 operator-(const vec3& other) const;
    vec3 operator-() const;
    vec3& operator+=(const vec3& other);
    bool operator==(const vec3& other) const;
    double operator*(const vec3& other) const;
    vec3 operator*(const int& other) const;
    vec3& operator/=(const double& other);
    float magnitude();
    vec3 update_by(vec3* v, double delta);
    vec3 direction(vec3 v3);
    vec3 normalize();
    double sqr_magn() const;
    void reset();
};

vec2 operator*(double scalar, vec2 v);

vec3 operator*(double scalar, vec3 v);
}