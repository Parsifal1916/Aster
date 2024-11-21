#ifndef PHYSICS_H
#define PHYSICS_H


double magnitude(vec2 v){
       return std::sqrt(v.x * v.x + v.y * v.y) ;
}

vec2 direction(vec2 v1, vec2 v2){
       return -v1 + v2;
}

vec2 normalize(vec2 v){
       double magn = magnitude(v);
       return vec2(v.x/magn, v.y/magn);
}

#endif