#include "Aster/impl/barnes_hut.tpp"
#include "Aster/impl/BHT_impl.tpp"

namespace Aster{
namespace Barnes{

double variation(){
    return ((std::rand() % 100) + 50)/150; 
}

template <>
vec3 pick_newpos(Simulation<vec3>* _s){
    return _s -> get_center() * variation();
}

template <>
vec2 pick_newpos(Simulation<vec2>* _s){
    return _s -> get_center() * variation();
}


template <>
size_t Barnes_Hut<vec2>::subdivide(int n){
    size_t child_index = nodes.size();
    Node<vec2>& cnode = nodes[n]; 

    cnode.child = child_index;

     for (int i = 0; i < num_childs; i++){
        nodes.push_back(Node<vec2>(
            {
            /*x = */ (0.5 - (i % 2 == 0)) * cnode.size + cnode.center.x,
            /*y = */ (0.5 - (i < 2)) * cnode.size + cnode.center.y
            },
            cnode.size/2
        ));
    }

    return child_index;
}

template <>
size_t Barnes_Hut<vec3>::subdivide(int n){
    size_t child_index = nodes.size();
    Node<vec3>& cnode = nodes[n]; 

    cnode.child = child_index;

    for (int i = 0; i < num_childs; i++){
       nodes.emplace_back(Node<vec3>(
           {
           /*x = */ (0.5 - (i % 2 == 0)) * cnode.size + cnode.center.x,
           /*y = */ (0.5 - (i < 2 || (i > 3 && i < 5))) * cnode.size + cnode.center.y,
           /*z = */ (0.5 - (i < 3)) * cnode.size + cnode.center.z 
           },
           cnode.size/2
       ));
    }

    return child_index;
}

/*
* returns the index of the child node
* given the position of the object to 
* insert and the center of the root node
*
* @param p: position of the object
* @param center: center of the node
* @retval: displacement in the node vector to find the optimal node
*/
template <>
int Barnes_Hut<vec2>::opt_position(vec2 p, vec2 center){ 
    bool 
        x = p.x > center.x,
        y = p.y > center.y
    ;

    return x + y*2;
}

/*
* returns the index of the child node
* given the position of the object to 
* insert and the center of the root node
*
* @param p: position of the object
* @param center: center of the node
* @retval: displacement in the node vector to find the optimal node
*/
template <>
int Barnes_Hut<vec3>::opt_position(vec3 p, vec3 center){ 
    bool 
        x = p.x > center.x,
        y = p.y > center.y,
        z = p.z > center.z
    ;

    return x + y*2 + z*4;
}

}
}