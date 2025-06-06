#pragma once

#include <cmath>

#include "Aster/physics/tool-chain.h"
#include "Aster/physics/vectors.h"

#include "Aster/simulations/sim_obj.h"
#include "Aster/simulations/basic.h"

#include "Aster/impl/config.h"

// NOTE: i don't trust loop unrolling

namespace Aster{ 

template <typename T> 
void update_SABA1(Simulation<T>* _s){
    update_euler(_s);
}

template <typename T> 
void update_SABA2(Simulation<T>* _s){
    constexpr double c1 = 0.21132486540518713447056597942719236;
    constexpr double c2 = 0.57735026918962573105886804114561528;
    constexpr double d1 = .5;
    const double dt = _s -> get_dt();

    for_each_body(_s, [c1, dt, _s](size_t body){
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c1 * dt;
    });

    _s -> update_forces(_s);

    for_each_body(_s, [d1, c2, dt, _s](size_t body){
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d1 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c2 * dt;        
    });

    _s -> update_forces(_s);

    for_each_body(_s, [d1, c1, dt, _s](size_t body){
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d1 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c1 * dt;
    });
}

template <typename T> 
void update_SABA3(Simulation<T>* _s){
    constexpr double c1 = 0.112701665379258297861042592558078468;
    constexpr double c2 = 0.387298334620741702138957407441921532;
    constexpr double d1 = 5.0/18.0;
    constexpr double d2 = 4.0/9.0;

    const double dt = _s -> get_dt();

    for_each_body(_s, [c1, dt, _s](size_t body){
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c1 * dt;
    });

    _s -> update_forces(_s);

    for_each_body(_s, [d1, c2, dt, _s](size_t body){
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d1 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c2 * dt;        
    });

    _s -> update_forces(_s);

    for_each_body(_s, [d2, c2, dt, _s](size_t body){
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d2 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c2 * dt;
    });

    _s -> update_forces(_s);

    for_each_body(_s, [d1, c1, dt, _s](size_t body){
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d1 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c1 * dt;
    });
}

template <typename T> 
void update_SABA4(Simulation<T>* _s){
    constexpr double c1 = 0.069431844202973713731097404888714663;
    constexpr double c2 = 0.260577634004598102102079337782924995;
    constexpr double c3 = 0.339981043584856257311344052141066641;

    constexpr double d1 = 0.173927422568726924856363780236279126;
    constexpr double d2 = 0.326072577431273047388060604134807363;

    const double dt = _s -> get_dt();

    for_each_body(_s, [c1, dt, _s](size_t body){
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c1 * dt;
    });

    _s -> update_forces(_s);

    for_each_body(_s, [d1, c2, dt, _s](size_t body){
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d1 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c2 * dt;        
    });

    _s -> update_forces(_s);

    for_each_body(_s, [d2, c3, dt, _s](size_t body){
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d2 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c3 * dt;        
    });

    _s -> update_forces(_s);

    for_each_body(_s, [d2, c2, dt, _s](size_t body){
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d2 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c2 * dt;        
    });

    _s -> update_forces(_s);

    for_each_body(_s, [d1, c1, dt, _s](size_t body){
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d1 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c1 * dt;        
    });
}

template <typename T>
void update_SABA5(Simulation<T>* _s) {
    constexpr double c1 = 0.046910077030668018149839326724759303;
    constexpr double c2 = 0.183855267916490455748501631205726881;
    constexpr double c3 = 0.269234655052841553857234657698427327;

    constexpr double d1 = 0.118463442528094542449679238416138105;
    constexpr double d2 = 0.239314335249683235451456653208879288;
    constexpr double d3 = 64.0 / 225.0;

    const double dt = _s -> get_dt();

    for_each_body(_s, [c1, dt, _s](size_t body) {
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c1 * dt;
    });

    _s -> update_forces(_s);

    for_each_body(_s, [d1, c2, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d1 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c2 * dt;
    });

    _s -> update_forces(_s);

    for_each_body(_s, [d2, c3, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d2 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c3 * dt;
    });

    _s -> update_forces(_s);

    for_each_body(_s, [d3, c3, dt, _s](size_t body) { 
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d3 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c3 * dt;
    });

    _s -> update_forces(_s);

    for_each_body(_s, [d2, c2, dt, _s](size_t body) { 
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d2 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c2 * dt;
    });

    _s -> update_forces(_s);

    for_each_body(_s, [d1, c1, dt, _s](size_t body) { 
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d1 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c1 * dt;
    });
}

template <typename T>
void update_SABA6(Simulation<T>* _s) {
    constexpr double c1 = 0.033765242898423986093849222753002695; 
    constexpr double c2 = 0.135630063868443757075450979737044631; 
    constexpr double c3 = 0.211295100191533802515448936669596706; 
    constexpr double c4 = 0.238619186083196908630501721680711935;

    constexpr double d1 = 0.085662246189585172520148071086366447;
    constexpr double d2 = 0.180380786524069303784916756918858056;
    constexpr double d3 = 0.233956967286345523694935171994775497;

    const double dt = _s -> get_dt();

    for_each_body(_s, [c1, dt, _s](size_t body) {
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c1 * dt;
    });

    _s -> update_forces(_s);

    for_each_body(_s, [d1, c2, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d1 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c2 * dt;
    });

    _s -> update_forces(_s);

    for_each_body(_s, [d2, c3, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d2 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c3 * dt;
    });

    _s -> update_forces(_s);

    for_each_body(_s, [d3, c4, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d3 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c4 * dt;
    });

    _s -> update_forces(_s);

    for_each_body(_s, [d2, c3, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d3 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c3 * dt;
    });

    _s -> update_forces(_s);

    for_each_body(_s, [d2, c2, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d2 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c2 * dt;
    });

    _s -> update_forces(_s);

    for_each_body(_s, [d1, c1, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d1 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c1 * dt;
    });

}

template <typename T>
void update_SABA7(Simulation<T>* _s) {
    constexpr double c1 = 0.025446043828620737736905157976074369; 
    constexpr double c2 = 0.103788363371682042331162455383531428; 
    constexpr double c3 = 0.167843017110998636478629180601913472; 
    constexpr double c4 = 0.202922575688698583453303206038480732; 

    constexpr double d1 = 0.064742483084434846635305716339541009;
    constexpr double d2 = 0.139852695744638333950733885711889791;  
    constexpr double d3 = 0.190915025252559472475184887744487567;
    constexpr double d4 = 256.0/1225.0;        

    const double dt = _s -> get_dt();

    for_each_body(_s, [c1, dt, _s](size_t body) {
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c1 * dt;
    });

    _s -> update_forces(_s);

    for_each_body(_s, [d1, c2, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d1 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c2 * dt;
    });

    _s -> update_forces(_s);

    for_each_body(_s, [d2, c3, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d2 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c3 * dt;
    });    

    _s -> update_forces(_s);

    for_each_body(_s, [d3, c4, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d3 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c4 * dt;
    });  

    _s -> update_forces(_s);

    for_each_body(_s, [d4, c3, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d4 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c4 * dt;
    });  

    _s -> update_forces(_s);

    for_each_body(_s, [d3, c2, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d3 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c3 * dt;
    });  

    _s -> update_forces(_s);

    for_each_body(_s, [d2, c1, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d2 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c2 * dt;
    }); 

    _s -> update_forces(_s);

    for_each_body(_s, [d2, c1, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d1 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c1 * dt;
    }); 
}
template <typename T>
void update_SABA8(Simulation<T>* _s) {
    constexpr double c1 = 0.019855071751231884158219565715263505; 
    constexpr double c2 = 0.081811689541954746046003466046821277; 
    constexpr double c3 = 0.135567033748648876886907443643292044; 
    constexpr double c4 = 0.171048883710339590439131453414531184;
    constexpr double c5 = 0.183434642495649804939476142360183981;

    constexpr double d1 = 0.050614268145188129576265677154981095;
    constexpr double d2 = 0.111190517226687235272177997213120442;
    constexpr double d3 = 0.156853322938943643668981100993300657;
    constexpr double d4 = 0.18134189168918099148257522463859781;

    const double dt = _s -> get_dt();

    for_each_body(_s, [c1, dt, _s](size_t body) {
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c1 * dt;
    });

    _s -> update_forces(_s);

    for_each_body(_s, [d1, c2, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d1 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c2 * dt;
    });

    _s -> update_forces(_s);

    for_each_body(_s, [d2, c3, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d2 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c3 * dt;
    });    

    _s -> update_forces(_s);

    for_each_body(_s, [d3, c4, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d3 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c4 * dt;
    });  

    _s -> update_forces(_s);

    for_each_body(_s, [d4, c5, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d4 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c5 * dt;
    });  

    _s -> update_forces(_s);

    for_each_body(_s, [d3, c4, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d4 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c4 * dt;
    });  

    _s -> update_forces(_s);

    for_each_body(_s, [d2, c3, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d3 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c3 * dt;
    });  

    _s -> update_forces(_s);

    for_each_body(_s, [d1, c2, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d2 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c2 * dt;
    });  

    _s -> update_forces(_s);

    for_each_body(_s, [d1, c2, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d1 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c1 * dt;
    });  
}

template <typename T>
void update_SABA9(Simulation<T>* _s) {
    constexpr double c1 = 0.015919880246186955082211898548163565; 
    constexpr double c2 = 0.066064566090495147768073207416968997; 
    constexpr double c3 = 0.111329837313022698495363874364130346; 
    constexpr double c4 = 0.144559004648390734135082012349068788; 
    constexpr double c5 = 0.162126711701904464519269007321668304;

    constexpr double d1 = 0.040637194180787205985946079055261825;
    constexpr double d2 = 0.090324080347428702029236015621456405;
    constexpr double d3 = 0.130305348201467731159371434709316425;
    constexpr double d4 = 0.156173538520001420034315203292221833;
    constexpr double d5 = 16384.0/99225.0;

    const double dt = _s -> get_dt();

    for_each_body(_s, [c1, dt, _s](size_t body) {
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c1 * dt;
    });

    _s -> update_forces(_s);

    for_each_body(_s, [d1, c2, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d1 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c2 * dt;
    });

    _s -> update_forces(_s);

    for_each_body(_s, [d2, c3, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d2 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c3 * dt;
    });    

    _s -> update_forces(_s);

    for_each_body(_s, [d3, c4, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d3 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c4 * dt;
    });  

    _s -> update_forces(_s);

    for_each_body(_s, [d4, c5, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d4 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c5 * dt;
    });   

    _s -> update_forces(_s);

    for_each_body(_s, [d5, c4, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d5 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c5 * dt;
    });   

    _s -> update_forces(_s);

    for_each_body(_s, [d4, c3, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d4 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c4 * dt;
    });  

    _s -> update_forces(_s);

    for_each_body(_s, [d3, c2, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d3 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c3 * dt;
    });  

    _s -> update_forces(_s);

    for_each_body(_s, [d2, c1, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d2 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c2 * dt;
    });  

    _s -> update_forces(_s);

    for_each_body(_s, [d2, c1, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d1 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c1 * dt;
    }); 
}


template <typename T>
void update_SABA10(Simulation<T>* _s) {
    constexpr double c1 = 0.013046735741414139961017993957773973; 
    constexpr double c2 = 0.054421580914093604672933661830479502; 
    constexpr double c3 = 0.092826899194980052248884661654309736; 
    constexpr double c4 = 0.123007087084888607717530710974544707; 
    constexpr double c5 = 0.142260527573807989957219971018032089; 
    constexpr double c6 = 0.148874338981631210884826001129719985;

    constexpr double d1 = 0.033335672154344068796784404946665896;
    constexpr double d2 = 0.074725674575290296572888169828848666;
    constexpr double d3 = 0.109543181257991021997767467114081596;
    constexpr double d4 = 0.134633359654998177545613460784734677;
    constexpr double d5 = 0.147762112357376435086946497325669165;

    const double dt = _s -> get_dt();

    for_each_body(_s, [c1, dt, _s](size_t body) {
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c1 * dt;
    });

    _s -> update_forces(_s);

    for_each_body(_s, [d1, c2, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d1 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c2 * dt;
    });

    _s -> update_forces(_s);

    for_each_body(_s, [d2, c3, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d2 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c3 * dt;
    });

    _s -> update_forces(_s);

    for_each_body(_s, [d3, c4, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d3 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c4 * dt;
    });

    _s -> update_forces(_s);

    for_each_body(_s, [d4, c5, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d4 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c5 * dt;
    });

    _s -> update_forces(_s);

    for_each_body(_s, [d5, c6, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d5 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c6 * dt;
    });

    _s -> update_forces(_s);

    for_each_body(_s, [d4, c5, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d5 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c5 * dt;
    });

    _s -> update_forces(_s);

    for_each_body(_s, [d3, c4, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d4 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c4 * dt;
    });

    _s -> update_forces(_s);

    for_each_body(_s, [d2, c3, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d3 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c3 * dt;
    });

    _s -> update_forces(_s);

    for_each_body(_s, [d1, c2, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d2 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c2 * dt;
    });

    _s -> update_forces(_s);

    for_each_body(_s, [d1, c2, dt, _s](size_t body) {
        _s -> bodies.get_velocity_of(body) += _s -> bodies.get_acc_of(body) * d1 * dt;
        _s -> bodies.get_position_of(body) += _s -> bodies.get_velocity_of(body) * c1 * dt;
    });
}
}