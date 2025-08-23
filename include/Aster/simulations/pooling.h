#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>


template <typename T>
struct memory_pool{
    size_t current_size;
    size_t used_memory = 0;

    void* data;

    public:
    void double_size();
    void reset();

    void* get_place(size_t size);

    inline T* operator[](size_t i){
        return (T*)this -> data + i;
    }

    ~memory_pool() {
        free(this->data);
    }
    template <typename... Args>
    T* insert(Args&&... args){
        auto* dict = (T*)this -> get_place(sizeof(T));

        return new (dict) T(std::forward<Args>(args)...);
    }

    inline size_t size(){
        return this -> used_memory;
    }

    memory_pool* make_memory_pool(size_t size) {
    memory_pool* pool = new memory_pool();

    pool -> sectors.reserve(10);
    pool -> sectors.push_back(0);
    pool -> data = malloc(size);
    pool -> current_size = size;
    pool -> used_memory = 0;

    return pool;
}


};

#include "pooling_helper.h"