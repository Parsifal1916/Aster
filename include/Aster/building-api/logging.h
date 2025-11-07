#pragma once 

#include <iostream>
#include <string>

namespace Aster{

enum error_type: int {
    LOW_t = 0,
    DEBUG_t = 1,
    WARNING_t = 2,
    ERROR_t = 3,
    CRITICAL_t = 4
};

constexpr int time_characters = 7;
extern error_type error_level;

double get_time();
bool critical_if(bool cond, std::string msg);
bool warn_if(bool cond, std::string msg);
bool err_if(bool cond, std::string msg);
void log_info(std::string msg, std::string EOL = "\n");

extern const std::string ERROR;
extern const std::string DEBUG;
extern const std::string CHECK;
extern const std::string WARNING;
extern const std::string CRITICAL;
extern const std::string RESET;

}