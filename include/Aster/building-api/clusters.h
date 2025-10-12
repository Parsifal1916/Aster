#pragma once 

#include <vector>
#include <mutex>
#include <string>
#include <functional>

#include "Aster/physics/vectors.h"
#include "Aster/physics/body.h"

namespace Aster{

class Simulation;
struct Cluster;
using cluster_builder_fptr = std::function<void(Cluster&, size_t)>;


struct Cluster{
    vec3 position; // cluster position 
    vec3 size; // cluster sidelenghts
    vec3 rotation; // rotation over all axis
    size_t number; // number of objects to load
    size_t loaded; // number of loaded objects
    REAL avr_mass; // average mass of the bodies

    public:
    // cluster name
    std::string name = "Custom Cluster 3d";

    // function to invoke to generate the bodies
    cluster_builder_fptr builder = nullptr;

    /**
    * @brief updates the number of loaded bodies by n
    * @param n: number of newly loaded bodies
    * @return a pointer to the cluster
    */
    Cluster* update(size_t n);


    /**
    * @brief returns the number of loaded bodies
    * @return the number of loaded bodies
    */
    size_t get_status();

    /**
    * @brief returns the number of bodies to load
    * @return the number of bodies to load
    */
    size_t get_objects();


    /**
    * @brief adds angles to the rotation parameter
    * @param x: rotation over x
    * @param y: rotation over y 
    * @param z: rotation over z
    * @return a pointer to the cluster
    */
    Cluster* rotate(REAL x, REAL y, REAL z = 0);

    /**
    * @brief adds angles to the rotation parameter
    * @param rot: rotation vector to add
    * @return a pointer to the cluster
    */
    Cluster* rotate(vec3 rot);


    /**
    * @brief moves the cluster
    * @param x: traslation over x
    * @param y: traslation over y 
    * @param z: traslation over z
    * @return a pointer to the cluster
    */
    Cluster* move(REAL x, REAL y, REAL z = 0);

    /**
    * @brief moves the cluster
    * @param dis: traslation vector
    * @return a pointer to the cluster
    */
    Cluster* move(vec3 dis);


    /**
    * @brief rotates the cluster
    * @param x: x rotation
    * @param y: y rotation 
    * @param z: z rotation
    * @return a pointer to the cluster
    */
    Cluster* set_rotation(REAL x, REAL y, REAL z = 0);

    /**
    * @brief rotates the cluester
    * @param rot: new rotation 
    * @return a pointer to the cluster
    */
    Cluster* set_rotation(vec3 rot); 

    /**
    * @brief moves the cluster
    * @param x: new x position 
    * @param y: new y position  
    * @param z: new z position 
    * @return a pointer to the cluster
    */
    Cluster* set_position(REAL x, REAL y, REAL z = 0);

    /**
    * @brief moves the cluester
    * @param pos: new position 
    * @return a pointer to the cluster
    */
    Cluster* set_position(vec3 pos);    
};

/** 
* @brief queue of clusters
* it is used to store unloaded clusters to then load them using built-in methods 
*/

struct ClusterQueue{
    // the actual vector of clusters
    std::vector<Cluster> data = {};

    public:
    std::mutex mtx;

    /**
    * @brief adds a cluster to the queue
    * @param cl: cluster to add
    * @return returns a pointer to the queue
    */
    ClusterQueue* add_cluster(Cluster cl);

    /**
    * @brief removes the last cluster
    * @return a pointer to the queue
    */
    ClusterQueue* pop_cluster();

    /**
    * @brief loads every cluster onto the given simulation
    * @param _s: simulation to load the bodies to
    * @return a pointer to the queue
    */
    ClusterQueue* load(Simulation* _s);

}; 
}