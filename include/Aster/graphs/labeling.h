#pragma once

#include <functional>
#include <string>

#include "Aster/simulations/sim_obj.h"

namespace Aster{
namespace Text{



using label_gen = std::function<std::string(Simulation*)>;
    
/**
* @brief struct used to generalize on screen labels
*/

struct Label{
    std::string content;
    label_gen generator;

    public:
    Label(label_gen gen) : generator(gen) {}
    
    /**
    * @brief updates the label's content using the generator function
    * @param _s: simulation to pass to the generator
    */
    void update(Simulation* _s);

    /**
    * @brief returns the content
    * @returns the content as a string
    */
    std::string get_content() const;
};


struct LabelQueue {
    // label queue
    std::vector<Label> queue;
    
    public:
    /**
    * @brief returns the lenght of the queue
    * @returns the lenght of the queue
    */
    size_t get_len() const;

    /**
    * @brief adds a label to the queue
    * @param gen: generator function for the labels (std::string(Simulation*))
    * @returns a pointer to the queue
    */
    LabelQueue* add_label(label_gen gen);

    /**
    * @brief clear the loading queue
    * @returns a pointer to the queue
    */
    LabelQueue* empty_queue();

    /**
    * @brief loads the labels onto the current window
    * @returns a pointer to the queue
    */
    void load(Simulation* _s);
};

}
}