#pragma once

#include "../sim-setup/simulation.h"
#include "vectors.h"
#include <SFML/Graphics.hpp>

using namespace Simulation;

/*
* gets the delta a from the previous 
* value of the scale factor 
*/
float get_a(float s){
       /* 
       * from friedmann we know that
       * a(t) = a_i exp(t * H_0 * root(omega))
       * we can differentiate and get
       * da = a_i * H_0 * root(omega) * exp ... *dt 
       */

       return initial_a *H_0 * root_omega * std::exp(time_passed* H_0 * root_omega);
}

/*
* updates the scale factor with the
* rk4 method
*/
void update_scale(){
       double k1 = get_a(current_a) * dt;
       double k2 = get_a(current_a + k1 * 0.5) * dt;
       double k3 = get_a(current_a + k2 * 0.5) * dt;
       double k4 = get_a(current_a + k3) * dt;

       current_a += (k1 + 2 * k2 + 2 * k3 + k4) / 6;
}

/*
* evaluates the post-newtonian approx.
* @param m1: mass of the first object
* @param m2: mass of the second object
* @param v1: velocity of the first object
* @param v2: velocity of the second object 
* @param p1: position of the first object
* @param p2: position if the second object
*/
vec2 pn2(double m1, double m2, vec2 v1, vec2 v2, vec2 p1, vec2 p2){
       
       vec2 n = ((p2 - p1) * current_a); 
       double r = (n.magnitude() + 1)*current_a;

       double term1 = G*m2/(c*c*r*r);
       double m_ratio = m1/m2;
       vec2 v = v1 - v2;

       double a_newton = G/(r*r + (double)e_squared);
       double correction = std::exp(-r*r/takeover)*divergence;
       vec2 pn1 = ( n*(4 + 2*m_ratio)*term1*c_squared*(double)r - v*v*n*(1+ 3*m_ratio) + v*n*v*4 - v*n*n*3/2) * term1 ;
       vec2 pn25 = (v*n)*v*-8/5 * G*G*m1*m2*(m1+m2)/(c_squared*c_squared*c*r*r*r);

       return n*( (double)a_newton) + pn1 + pn25;
}
std::vector<vec2> get_new_pos(vec2* position, vec2* velocity, vec2* acceleration, double step){
       std::vector<vec2> retval;
       retval.push_back(vec2({velocity -> x + acceleration -> x *step, velocity -> y + acceleration -> y *step}));
       retval.push_back(vec2({position -> x + retval[0].x * step     , position -> y + retval[0].y * step     }));
       return retval;
}
vec2 rk4(double m1, double m2, vec2 v1, vec2 v2, vec2 p1, vec2 p2){
       vec2 k1, k2, k3, k4;
       std::vector<vec2> dummy;
       k1 = pn2(m1, m2, v1, v2, p1, p2) * dt;
       dummy = get_new_pos(&p1, &v1, &k1, dt/2);
       k2 = pn2(m1, m2, v1 + dummy[0]/2.f, v2, p1 + dummy[1]/2.f, p2) * dt;
       dummy = get_new_pos(&p1, &v1, &k2, dt/2);
       k3 = pn2(m1, m2, v1 + dummy[0]/2.f, v2, p1 + dummy[1]/2.f, p2) * dt;
       dummy = get_new_pos(&p1, &v1, &k3, dt);
       k4 = pn2(m1, m2, v1 + dummy[0]/2.f, v2, p1 + dummy[1]/2.f, p2) * dt;

       return (k1 + k2*2.f + k3*2.f + k4)/6.f;
}



struct Body{
public:
       double mass = 1;
       double temp = 0;
       vec2 position = {};
       vec2 velocity = {};
       vec2 acceleration = {};
       sf::Color color = sf::Color::White; 
       unsigned int index = 0;

       Body(double mass, vec2 position, vec2 velocity, unsigned int index, sf::Color c ) : color(c),mass(mass), index(index), position(position), velocity(velocity){};
       Body(){}

       static sf::Color get_color(double t) {
              uint8_t r, g, b;

              if (t <= 0.3) {
                     double normalizedT = t / 0.33;
                     r = 0 + normalizedT * (255 - 0);
                     g = 0;
                     b = 255 - normalizedT * (255 - 0);
              } else{
                     double normalizedT = (t - 0.33) / 0.33;
                     r = 255;
                     g = 0 + normalizedT * (255 - 0);
                     b = 0;
              }
              return {r, g, b};
       }

       static double sigmoid(double x) {
              return 1.0 / (1.0 + exp(-std::log10(x)));
       }
};
