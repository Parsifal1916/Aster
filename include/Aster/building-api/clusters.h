#pragma once 

#include <vector>
#include <mutex>
#include <string>
#include <functional>

#include "Aster/physics/vectors.h"
#include "Aster/physics/body.h"


namespace Aster{
template <typename T> struct Cluster;
template <typename T> struct Simulation;

template <typename T> 
using cluster_builder_fptr = std::function<Body<T>(Cluster<T>&, size_t)>;

template <typename T> 
struct Cluster{
    T position;
    T size;
    T rotation;
    size_t number;
    size_t loaded;
    double avr_mass;

    public:
    
    std::string name = "Custom Cluster 3d";

    cluster_builder_fptr<T> builder = nullptr;

    Cluster<T>* update(size_t n);

    size_t get_status();
    size_t get_objects();

    Cluster<T>* rotate(double x, double y, double z = 0);
    Cluster<T>* rotate(T rot);

    Cluster<T>* move(double x, double y, double z = 0);
    Cluster<T>* move(T dis);

    Cluster<T>* set_rotation(double x, double y, double z = 0);
    Cluster<T>* set_rotation(T rot); 

    Cluster<T>* set_position(double x, double y, double z = 0);
    Cluster<T>* set_position(T pos);    
};

template <typename T>
struct ClusterQueue{
    std::vector<Cluster<T>> data = {};

    public:
    std::mutex mtx;

    ClusterQueue<T>* add_cluster(Cluster<T> cl2d);
    ClusterQueue<T>* pop_cluster();
    ClusterQueue<T>* load(Simulation<T>* _s);

}; 
}

#include "Aster/impl/clusters_endpoint.tpp"