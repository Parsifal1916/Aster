#include "Aster/graphs/labeling.h"
#include "Aster/graphics/text_rendering.h"
#include "Aster/building-api/logging.h"

namespace Aster{
namespace Text{

/**
* @brief updates the label's content using the generator function
*/

void Label::update(Simulation* _s){
    content = generator(_s); // executes the function onto the simualtion
    warn_if(content == "", "no text generated for label");
}

/**
* @brief returns the content
* @returns the content as a string
*/

std::string Label::get_content() const{
    return content;
}

/**
* @brief returns the lenght of the queue
* @returns the lenght of the queue
*/

size_t LabelQueue::get_len() const{
    return this -> queue.size();
};

/**
* @brief adds a label to the queue
* @param gen: generator function for the labels (std::string(Simulation*))
* @returns a pointer to the queue
*/

LabelQueue* LabelQueue::add_label(label_gen gen){
    if (warn_if(get_len() > 10, "cannot add label, maximum capacity of 10 reached"))
        return this;

    this -> queue.push_back(Label(gen));
    return this;
};

/**
* @brief clear the loading queue
* @returns a pointer to the queue
*/

LabelQueue* LabelQueue::empty_queue(){
    if (warn_if(!get_len(), "cannot empty label queue, it already is empty!"))
        return this;

    bool _ = this -> queue.empty();
    return this;
}

/**
* @brief loads the labels onto the current window
* @returns a pointer to the queue
*/

void LabelQueue::load(Simulation* _s){
    constexpr int
        font_size = 15,
        spacing = 30,
        x_margin = 50,
        y_margin = 50
    ;

    for (int i = 0; i < get_len(); ++i){
        queue[i].update(_s);
        write2screen(x_margin, y_margin + i * spacing, font_size,queue[i].get_content().c_str());
    }

}

}
}