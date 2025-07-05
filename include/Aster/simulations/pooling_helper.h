#include <cassert>
#include <new>
/*
* genera la pool e ne ritorna un pointer
* @param size: grandezza della pool
* @return: un ptr alla pool
*/
template <typename T>
memory_pool<T>* make_memory_pool(size_t size) {
    memory_pool<T>* pool = new memory_pool<T>();

    pool -> data = malloc(size);
    pool -> current_size = size;
    pool -> used_memory = 0;

    return pool;
}

template <typename T>
void memory_pool<T>::REAL_size() {
    current_size *= 2;
    data = std::aligned_alloc(alignof(T), current_size); // lento
}

template <typename T>
void memory_pool<T>::reset(){
    used_memory = 0;
}

template <typename T>
void* memory_pool<T>::get_place(size_t size) {
    if (this -> used_memory + size > this -> current_size) 
        this -> REAL_size();

    void* ptr = (char*)this -> data + this -> used_memory;
    this -> used_memory += size;
    return ptr;
}
