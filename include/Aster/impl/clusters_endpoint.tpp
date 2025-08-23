#pragma once

#include "Aster/building-api/clusters.h"

#include "Aster/physics/vectors.h"
#include "Aster/physics/body.h"
#include "Aster/building-api/logging.h"

namespace Aster{

//===---------------------------------------------------------===//
// 3d clusters impl                                              //
//===---------------------------------------------------------===//

/**
* @brief updates the number of loaded bodies by n
* @param n: number of newly loaded bodies
* @return a pointer to the cluster
*/
template <typename T>
Cluster<T>* Cluster<T>::update(size_t n){
    loaded += n;
    return this;
} 

/**
* @brief returns the number of loaded bodies
* @return the number of loaded bodies
*/
template <typename T>
size_t Cluster<T>::get_status(){
    return loaded;
} 

/**
* @brief returns the number of bodies to load
* @return the number of bodies to load
*/
template <typename T>
size_t Cluster<T>::get_objects(){
    return number;
} 

/**
* @brief changes the rotation parameter
* @param x: rotation over x
* @param y: rotation over y   
* @param z: rotation over z
* @return a pointer to the cluster
*/
template <typename T>
Cluster<T>* Cluster<T>::rotate(REAL x, REAL y, REAL z){
    rotation += vec3({x, y, z});
    return this;
} 

/**
* @brief adds angles to the rotation parameter
* @param rot: rotation vector to add
* @return a pointer to the cluster
*/
template <typename T>
Cluster<T>* Cluster<T>::rotate(T rot){
    rotation += rot;
    return this;
} 


/**
* @brief moves the cluster
* @param x: traslation over x
* @param y: traslation over y 
* @param z: traslation over z
* @return a pointer to the cluster
*/
template <typename T>
Cluster<T>* Cluster<T>::move(REAL x, REAL y, REAL z){
    position += vec3({x, y, z});
    return this;
} 

/**
* @brief moves the cluster
* @param dis: traslation vector
* @return a pointer to the cluster
*/
template <typename T>
Cluster<T>* Cluster<T>::move(T dis){
    position += dis;
    return this;
} 

/**
* @brief rotates the cluster
* @param x: x rotation
* @param y: y rotation 
* @param z: z rotation
* @return a pointer to the cluster
*/
template <typename T>
Cluster<T>* Cluster<T>::set_rotation(REAL x, REAL y, REAL z){
    rotation = vec3({x, y, z});
    return this;
} 

/**
* @brief rotates the cluester
* @param rot: new rotation 
* @return a pointer to the cluster
*/
template <typename T>
Cluster<T>* Cluster<T>::set_rotation(T rot){
    rotation = rot;
    return this;
} 

/**
* @brief moves the cluster
* @param x: new x position 
* @param y: new y position  
* @param z: new z position 
* @return a pointer to the cluster
*/
template <typename T>
Cluster<T>* Cluster<T>::set_position(REAL x, REAL y, REAL z){
    position = vec3({x, y, z});
    return this;
} 

/**
* @brief moves the cluester
* @param  pos: new position 
* @return a pointer to the cluster
*/
template <typename T>
Cluster<T>* Cluster<T>::set_position(T pos){
    position = pos;
    return this;
} 


//===---------------------------------------------------------===//
// 2d queues impl                                                //
//===---------------------------------------------------------===//

/**
* @brief adds a cluster to the queue
* @param cl: cluster to add
* @return returns a pointer to the queue
*/
template <typename T>
ClusterQueue<T>* ClusterQueue<T>::add_cluster(Cluster<T> cl){
    // acquires the general lock before adding the cluster
    std::lock_guard<std::mutex> lock(mtx);

    data.emplace_back(cl);
    return this;
}

/**
* @brief removes the last cluster
* @return a pointer to the queue
*/
template <typename T>
ClusterQueue<T>* ClusterQueue<T>::pop_cluster(){
    // checks if there are still elements in the queue
    if (warn_if(!data.size(), "queue has no data inside, cannot remove elements from it"))
        return this;

    // acquires the general lock before adding the cluster
    std::lock_guard<std::mutex> lock(mtx);

    // actually removes the last element and quits
    data.pop_back();
    return this;
}

/**
* @brief loads every cluster onto the given simulation
* @param _s: simulation to load the bodies to
* @return a pointer to the queue
*/
template <typename T>
ClusterQueue<T>* ClusterQueue<T>::load(Simulation<T>* _s){
    warn_if(data.size() == 0, "no clusters to load from the queue"); 

    if (critical_if(!_s, "simulation given is nullptr"))
        exit(-1);

    while (data.size()){ // takes the size as a condition since it is gradually shortening it
        // saves the last cluster to a reference
        auto& cluster = data.back();
        
        // gives feedback on what it is loading
        _s -> loading_meta.first = cluster.name;
        for (size_t i = 0; i < cluster.number; ++i){
            cluster.builder(cluster, cluster.get_status());

            // updates the number of loaded bodies and the percentage
            cluster.update(1);
            _s -> loading_meta.second = cluster.get_status() / cluster.number;

        }

        data.pop_back(); // pops the last element    
    }

    return this;
}


}
