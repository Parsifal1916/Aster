#pragma once

#include "Aster/building-api/clusters.h"

#include "Aster/physics/vectors.h"
#include "Aster/physics/body.h"

namespace Aster{

//===---------------------------------------------------------===//
// 3d clusters impl                                              //
//===---------------------------------------------------------===//

template <typename T>
size_t Cluster<T>::get_status(){
    return loaded;
} 

template <typename T>
size_t Cluster<T>::get_objects(){
    return number;
} 


template <typename T>
Cluster<T>* Cluster<T>::rotate(double x, double y, double z){
    rotation += vec3({x, y, z});
    return this;
} 

template <typename T>
Cluster<T>* Cluster<T>::rotate(T rot){
    rotation += rot;
    return this;
} 

template <typename T>
Cluster<T>* Cluster<T>::move(double x, double y, double z){
    position += vec3({x, y, z});
    return this;
} 

template <typename T>
Cluster<T>* Cluster<T>::move(T dis){
    position += dis;
    return this;
} 

template <typename T>
Cluster<T>* Cluster<T>::set_rotation(double x, double y, double z){
    rotation = vec3({x, y, z});
    return this;
} 

template <typename T>
Cluster<T>* Cluster<T>::set_rotation(T rot){
    rotation = rot;
    return this;
} 

template <typename T>
Cluster<T>* Cluster<T>::set_position(double x, double y, double z){
    position = vec3({x, y, z});
    return this;
} 

template <typename T>
Cluster<T>* Cluster<T>::set_position(T pos){
    position = pos;
    return this;
} 

template <typename T>
Cluster<T>* Cluster<T>::update(size_t n){
    loaded += n;
    return this;
} 

//===---------------------------------------------------------===//
// 2d queues impl                                                //
//===---------------------------------------------------------===//

template <typename T>
ClusterQueue<T>* ClusterQueue<T>::add_cluster(Cluster<T> cl){
    std::lock_guard<std::mutex> lock(mtx);

    data.emplace_back(cl);
    return this;
}

template <typename T>
ClusterQueue<T>* ClusterQueue<T>::pop_cluster(){
    std::lock_guard<std::mutex> lock(mtx);

    data.pop_back();
    return this;
}

template <typename T>
ClusterQueue<T>* ClusterQueue<T>::load(Simulation<T>* _s){
    while (data.size()){
        auto& cluster = data.back();
        
        _s -> loading_meta.first = cluster.name;
        for (size_t i = 0; i < cluster.number; ++i){
            _s -> bodies.emplace_back(
                cluster.builder(cluster, cluster.get_status())
            );
            cluster.update(1);
            _s -> loading_meta.second = cluster.get_status() / cluster.number;

        }

        data.pop_back();    
    }

    return this;
}


}
