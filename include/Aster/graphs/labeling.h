#pragma once

#include <functional>
#include <string>

#include "Aster/simulations/sim_obj.h"

namespace Aster{
namespace Text{


template <typename T>
using label_gen = std::function<std::string(Simulation<T>*)>;
    
/**
* @brief struct used to generalize on screen labels
*/
template<typename T>
struct Label{
    std::string content;
    label_gen<T> generator;

    public:
    Label(label_gen<T> gen) : generator(gen) {}
    
    /**
    * @brief updates the label's content using the generator function
    * @param _s: simulation to pass to the generator
    */
    void update(Simulation<T>* _s);

    /**
    * @brief returns the content
    * @returns the content as a string
    */
    std::string get_content() const;
};

template <typename T>
struct LabelQueue {
    // label queue
    std::vector<Label<T>> queue;
    
    public:
    /**
    * @brief returns the lenght of the queue
    * @returns the lenght of the queue
    */
    size_t get_len() const;

    /**
    * @brief adds a label to the queue
    * @param gen: generator function for the labels (std::string(Simulation<T>*))
    * @returns a pointer to the queue
    */
    LabelQueue<T>* add_label(label_gen<T> gen);

    /**
    * @brief clear the loading queue
    * @returns a pointer to the queue
    */
    LabelQueue<T>* empty_queue();

    /**
    * @brief loads the labels onto the current window
    * @returns a pointer to the queue
    */
    void load(Simulation<T>* _s);
};

}
}

#include "Aster/impl/labeling.tpp"