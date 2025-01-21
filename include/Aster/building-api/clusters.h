#pragma once 

#include <vector>
#include <mutex>
#include <string>
#include <functional>

#include "Aster/physics/vectors.h"
#include "Aster/physics/body.h"


namespace Aster{

using cluster_builder_fptr2d = std::function<Body(struct Cluster2d&, size_t)>;
using cluster_builder_fptr3d = std::function<Body3d(struct Cluster3d&, size_t)>;

struct Cluster3d{
    vec3 position = {0,0,0};
    vec3 size = {0,0,0};
    vec3 rotation = {0,0,0};
    size_t number;
    size_t loaded;
    double avr_mass;

    public:
    
    std::string name = "Custom Cluster 3d";

    cluster_builder_fptr3d builder = nullptr;

    Cluster3d* update(size_t n);

    size_t get_status();
    size_t get_objects();

    Cluster3d* rotate(double x, double y, double z);
    Cluster3d* rotate(vec3 rot);

    Cluster3d* move(double x, double y, double z);
    Cluster3d* move(vec3 dis);

    Cluster3d* set_rotation(double x, double y, double z);
    Cluster3d* set_rotation(vec3 rot); 

    Cluster3d* set_position(double x, double y, double z);
    Cluster3d* set_position(vec3 pos);    
};

struct Cluster2d{
    vec2 position = {0,0};
    vec2 size = {0,0};
    vec2 rotation = {0,0};
    size_t number;
    double avr_mass;

    size_t loaded = 0; 
    cluster_builder_fptr2d builder = nullptr;

    public:
    
    std::string name = "Custom Cluster 2d";

    Cluster2d* update(size_t n);

    size_t get_status();
    size_t get_objects();

    Cluster2d* move(double x, double y);
    Cluster2d* move(vec2 dis);

    Cluster2d* set_position(double x, double y);
    Cluster2d* set_position(vec2 rot); 

};

struct Queue2d{
    std::vector<Cluster2d> data = {};

    public:
    std::mutex mtx;

    Queue2d* add_cluster(Cluster2d cl2d);
    Queue2d* pop_cluster();
    Queue2d* load(struct Simulation* _s);

}; 

struct Queue3d{
    std::vector<Cluster3d> data = {};

    public:
    std::mutex mtx;

    Queue3d* add_cluster(Cluster3d cl3d);
    Queue3d* pop_cluster();
    Queue3d* load(struct Simulation3d* _s);

}; 


}