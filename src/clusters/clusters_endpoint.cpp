#include "Aster/building-api/clusters.h"

#include "Aster/physics/vectors.h"
#include "Aster/physics/body.h"

#include "Aster/simulations/sim_obj.h"
#include "Aster/simulations/3d_sim_obj.h"

namespace Aster{

//===---------------------------------------------------------===//
// 3d clusters impl                                              //
//===---------------------------------------------------------===//

size_t Cluster3d::get_status(){
    return loaded;
} 

size_t Cluster3d::get_objects(){
    return number;
} 


Cluster3d* Cluster3d::rotate(double x, double y, double z){
    rotation += vec3({x, y, z});
    return this;
} 

Cluster3d* Cluster3d::rotate(vec3 rot){
    rotation += rot;
    return this;
} 

Cluster3d* Cluster3d::move(double x, double y, double z){
    position += vec3({x, y, z});
    return this;
} 

Cluster3d* Cluster3d::move(vec3 dis){
    position += dis;
    return this;
} 

Cluster3d* Cluster3d::set_rotation(double x, double y, double z){
    rotation = vec3({x, y, z});
    return this;
} 

Cluster3d* Cluster3d::set_rotation(vec3 rot){
    rotation = rot;
    return this;
} 

Cluster3d* Cluster3d::set_position(double x, double y, double z){
    position = vec3({x, y, z});
    return this;
} 

Cluster3d* Cluster3d::set_position(vec3 pos){
    position = pos;
    return this;
} 

Cluster3d* Cluster3d::update(size_t n){
    loaded += n;
    return this;
} 


//===---------------------------------------------------------===//
// 2d clusters impl                                              //
//===---------------------------------------------------------===//

size_t Cluster2d::get_status(){
    return loaded;
} 

size_t Cluster2d::get_objects(){
    return number;
} 


Cluster2d* Cluster2d::move(double x, double y){
    position += vec2({x, y});
    return this;
} 

Cluster2d* Cluster2d::move(vec2 dis){
    position += dis;
    return this;
} 

Cluster2d* Cluster2d::set_position(double x, double y){
    position = vec2({x, y});
    return this;
} 

Cluster2d* Cluster2d::set_position(vec2 pos){
    position = pos;
    return this;
} 

Cluster2d* Cluster2d::update(size_t n){
    loaded += n;
    return this;
} 

//===---------------------------------------------------------===//
// 2d queues impl                                                //
//===---------------------------------------------------------===//

Queue2d* Queue2d::add_cluster(Cluster2d cl2d){
    std::lock_guard<std::mutex> lock(mtx);

    data.emplace_back(cl2d);
    return this;
}

Queue2d* Queue2d::pop_cluster(){
    std::lock_guard<std::mutex> lock(mtx);

    data.pop_back();
    return this;
}

Queue2d* Queue2d::load(Simulation* _s){
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

//===---------------------------------------------------------===//
// 2d queues impl                                                //
//===---------------------------------------------------------===//

Queue3d* Queue3d::add_cluster(Cluster3d cl3d){
    std::lock_guard<std::mutex> lock(mtx);

    data.emplace_back(cl3d);
    return this;
}

Queue3d* Queue3d::pop_cluster(){
    std::lock_guard<std::mutex> lock(mtx);

    data.pop_back();
    return this;
}

Queue3d* Queue3d::load(Simulation3d* _s){
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
