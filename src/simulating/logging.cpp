#include <iostream>
#include <chrono>
#include <iomanip>

#include "Aster/building-api/logging.h"

namespace Aster{

const std::string ERROR = "\033[31m";
const std::string DEBUG = "\033[34m";
const std::string CHECK = "\033[32m";
const std::string WARNING = "\033[33m";
const std::string CRITICAL = "\033[35m";
const std::string RESET = "\033[0m";


error_type error_level = WARNING_t;

double get_time(){
    using namespace std::chrono;
    return duration<double>(steady_clock::now().time_since_epoch()).count() / 1000000;
}

bool critical_if(bool cond, std::string msg){
    if (!cond) return false ;

    std::cout << "("
              << std::fixed
              << std::setprecision(3)
              << std::setw(time_characters) 
              << get_time() << "s"
              << ")" << "[" << CRITICAL << "CRIT" << RESET << "] " << msg << "\n";
    return true;
}

bool warn_if(bool cond, std::string msg){
    if (!cond) return false;

    std::cout << "("
              << std::fixed
              << std::setprecision(3)
              << std::setw(time_characters) 
              << get_time() << "s" 
              << ")" << "[" << WARNING << "WARN" << RESET << "] " << msg << "\n";
    return true;
}


bool err_if(bool cond, std::string msg){
    if (!cond) return false;

    std::cout << "("
              << std::fixed
              << std::setprecision(3)
              << std::setw(time_characters) 
              << get_time() << "s"
              << ")" << "[" << ERROR << "ERR" << RESET << "] " << msg << "\n";
    return true;
}


void log_info(std::string msg){

    std::cout << "("
              << std::fixed
              << std::setprecision(3)
              << std::setw(time_characters) 
              << get_time() << "s" 
              << ")" << "[" << DEBUG << "INFO" << RESET << "] " <<  msg << "\n";
}

}