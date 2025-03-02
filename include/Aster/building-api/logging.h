#pragma once 

#include <iostream>
#include <string>

namespace Aster{

enum error_type: int {
    DEBUG_t = 0,
    WARNING_t = 1,
    ERROR_t = 2,
    CRITICAL_t = 3
};

constexpr int time_characters = 7;
extern error_type error_level;

double get_time();
void critical_if(bool cond, std::string msg);
void warn_if(bool cond, std::string msg);
void err_if(bool cond, std::string msg);
void log_info(std::string msg);

extern const std::string ERROR;
extern const std::string DEBUG;
extern const std::string CHECK;
extern const std::string WARNING;
extern const std::string CRITICAL;
extern const std::string RESET;

}