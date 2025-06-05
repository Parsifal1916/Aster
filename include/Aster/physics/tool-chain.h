#pragma once
#include <vector>
#include <map>
#include <string>
#define CL_TARGET_OPENCL_VERSION 300
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS


#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif

#include "Aster/physics/body.h"
#include "Aster/simulations/basic.h"

namespace Aster{

template <typename T> class Simulation;

extern double c1;
extern double c2;
extern double d1;
extern double d2;

//===---------------------------------------------------------===//
// 2d methods                                                    //
//===---------------------------------------------------------===//

template <typename T> class Simulation;

template <typename T>
void update_euler(Simulation<T>* _s);

template <typename T>
void update_symplectic4(Simulation<T>* _s);

template <typename T>
void update_SABA1(Simulation<T>* _s);

template <typename T>
void update_SABA2(Simulation<T>* _s);

template <typename T>
void update_SABA3(Simulation<T>* _s);

template <typename T>
void update_SABA4(Simulation<T>* _s);

template <typename T>
void update_SABA4(Simulation<T>* _s);

template <typename T>
void update_SABA5(Simulation<T>* _s);

template <typename T>
void update_SABA6(Simulation<T>* _s);

template <typename T>
void update_SABA7(Simulation<T>* _s);

template <typename T>
void update_SABA8(Simulation<T>* _s);

template <typename T>
void update_SABA9(Simulation<T>* _s);

template <typename T>
void update_SABA10(Simulation<T>* _s);

template <typename T>
T newtonian(double m1, double m2, T v1, T v2, T p1, T p2, Simulation<T>* _s);

/*
* evaluates the post-newtonian approx.
* @param m1: mass of the first object
* @param m2: mass of the second object
* @param v1: velocity of the first object
* @param v2: velocity of the second object 
* @param p1: position of the first object
* @param p2: position if the second object
*/
template <typename T>
T pn2(double m1, double m2, T v1, T v2, T p1, T p2, Simulation<T>* _s);

template <typename T>
std::vector<T> get_new_pos(T* position, T* velocity, T* acceleration, double step);

template <typename T>
T rk4(double m1, double m2, T v1, T v2, T p1, T p2, Simulation<T>* _s);

template <typename T, typename F>
void for_each_body(Simulation<T>* _s, F func);

template <typename T>
double get_eccentricity(Simulation<T>* _s, size_t body, double relv_sq, double w_squared,  double radius, double mass2);

template <typename T>
void get_new_temp(Simulation<T>* _s, size_t body, size_t body2);

template <typename T>
void get_new_temp(Simulation<T>* _s, size_t body, T pos, T vel, double temp, double mass);

template <typename T>
void compute_rad_pressure(Simulation<T>* _s, size_t body, T pos, double temp);
 

//===---------------------------------------------------------===//
// Cosmology stuff                                               //
//===---------------------------------------------------------===//

/*
* gets the delta a from the previous 
* value of the scale factor 
*/
template <typename T>
float get_da(float s, Simulation<T>* _s);

/*
* updates the scale factor with the
* rk4 method
*/
template <typename T>
void update_scale(Simulation<T>* _s);

/*
*/

namespace GPU{
//===---------------------------------------------------------===//
// GPU optimized stuff                                           //
//===---------------------------------------------------------===//

cl_platform_id platform;
cl_device_id device;

cl_context context;
cl_command_queue queue;

bool has_initialized = false;

inline constexpr double SABA1_coeff[] = {0,
    2, 2, 2
};

inline constexpr double SABA2_coeff[] = {0,
    0.21132486540518713447056597942719236,
    .5, 
    0.57735026918962573105886804114561528, 
    .5, 
    0.21132486540518713447056597942719236
};

inline constexpr double SABA3_coeff[] = {0,
    0.112701665379258297861042592558078468,
    5.0/18.0,
    0.387298334620741702138957407441921532,
    4.0/9.0,
    0.387298334620741702138957407441921532,
    5.0/18.0,
    0.112701665379258297861042592558078468
}; 

inline constexpr double SABA4_coeff[] = {0,
    0.069431844202973713731097404888714663,
    0.173927422568726924856363780236279126,
    0.260577634004598102102079337782924995,
    0.326072577431273047388060604134807363,
    0.339981043584856257311344052141066641,
    0.326072577431273047388060604134807363,
    0.260577634004598102102079337782924995,
    0.173927422568726924856363780236279126,
    0.069431844202973713731097404888714663
}; 

inline constexpr double SABA5_coeff[] = {0,
    0.046910077030668018149839326724759303,
    0.118463442528094542449679238416138105,
    0.183855267916490455748501631205726881,
    0.239314335249683235451456653208879288,
    0.269234655052841553857234657698427327,
    64.0 / 225.0,
    0.269234655052841553857234657698427327,
    0.239314335249683235451456653208879288,
    0.183855267916490455748501631205726881,
    0.118463442528094542449679238416138105,
    0.046910077030668018149839326724759303
}; 

inline constexpr double SABA6_coeff[] = { 0,   
    0.033765242898423986093849222753002695, 
    0.085662246189585172520148071086366447,
    0.135630063868443757075450979737044631,
    0.180380786524069303784916756918858056,
    0.211295100191533802515448936669596706, 
    0.233956967286345523694935171994775497,
    0.238619186083196908630501721680711935,
    0.233956967286345523694935171994775497,
    0.211295100191533802515448936669596706,
    0.180380786524069303784916756918858056,
    0.135630063868443757075450979737044631,
    0.085662246189585172520148071086366447,
    0.033765242898423986093849222753002695
}; 

inline constexpr double SABA7_coeff[] = {0,
    0.025446043828620737736905157976074369, 
    0.064742483084434846635305716339541009,
    0.103788363371682042331162455383531428, 
    0.139852695744638333950733885711889791,
    0.167843017110998636478629180601913472, 
    0.190915025252559472475184887744487567,
    0.202922575688698583453303206038480732,
    256.0/1225.0,
    0.202922575688698583453303206038480732,
    0.190915025252559472475184887744487567,
    0.167843017110998636478629180601913472,
    0.139852695744638333950733885711889791,
    0.103788363371682042331162455383531428,
    0.064742483084434846635305716339541009,
    0.025446043828620737736905157976074369
}; 

inline constexpr double SABA8_coeff[] = {0,
    0.019855071751231884158219565715263505,  
    0.050614268145188129576265677154981095,
    0.081811689541954746046003466046821277, 
    0.111190517226687235272177997213120442,
    0.135567033748648876886907443643292044, 
    0.156853322938943643668981100993300657,
    0.171048883710339590439131453414531184,
    0.18134189168918099148257522463859781,    
    0.183434642495649804939476142360183981,
    0.18134189168918099148257522463859781,    
    0.171048883710339590439131453414531184,
    0.156853322938943643668981100993300657,
    0.135567033748648876886907443643292044,
    0.111190517226687235272177997213120442,
    0.081811689541954746046003466046821277,
    0.050614268145188129576265677154981095,
    0.019855071751231884158219565715263505
}; 

inline constexpr double SABA9_coeff[] = {0,
    0.015919880246186955082211898548163565, 
    0.040637194180787205985946079055261825,
    0.066064566090495147768073207416968997, 
    0.090324080347428702029236015621456405,
    0.111329837313022698495363874364130346, 
    0.130305348201467731159371434709316425,
    0.144559004648390734135082012349068788, 
    0.156173538520001420034315203292221833,
    0.162126711701904464519269007321668304,
    16384.0/99225.0,
    0.162126711701904464519269007321668304,
    0.156173538520001420034315203292221833,
    0.144559004648390734135082012349068788,
    0.130305348201467731159371434709316425,
    0.111329837313022698495363874364130346,
    0.090324080347428702029236015621456405,
    0.066064566090495147768073207416968997,
    0.040637194180787205985946079055261825,
    0.015919880246186955082211898548163565
}; 

inline constexpr double SABA10_coeff[] = { 0,   
    0.013046735741414139961017993957773973, 
    0.033335672154344068796784404946665896,
    0.054421580914093604672933661830479502, 
    0.074725674575290296572888169828848666,
    0.092826899194980052248884661654309736, 
    0.109543181257991021997767467114081596,
    0.123007087084888607717530710974544707, 
    0.134633359654998177545613460784734677,
    0.142260527573807989957219971018032089, 
    0.147762112357376435086946497325669165,
    0.148874338981631210884826001129719985,
    0.147762112357376435086946497325669165,
    0.142260527573807989957219971018032089,
    0.134633359654998177545613460784734677,
    0.123007087084888607717530710974544707,
    0.109543181257991021997767467114081596,
    0.092826899194980052248884661654309736,
    0.074725674575290296572888169828848666,
    0.054421580914093604672933661830479502,
    0.033335672154344068796784404946665896,
    0.013046735741414139961017993957773973
}; 

inline constexpr size_t saba_coeff_lng[10] = {4, 6, 8, 10, 12, 14, 16, 18, 20, 22};

inline constexpr const double* saba_coeffs[] = {
    SABA1_coeff, SABA2_coeff, 
    SABA3_coeff, SABA4_coeff, 
    SABA5_coeff, SABA6_coeff, 
    SABA7_coeff, SABA8_coeff,
    SABA9_coeff, SABA10_coeff 
};

/**
* @brief sets the best device in terms of compute power
*/
void select_best_device();

/**
* @brief compiles a kernel
* @param name: pointer to the name of the function
* @param source: source code of the kernel
* @param k: kernel object onto which to write the kernel 
*/
void compile_kernel(std::string* name, std::string* source, cl_kernel& k);

/**
* @brief initializes opencl and finds the right device
*/
void init_opencl();

/**
* @brief compiles the force program
*/
template <typename T>
func_ptr<T> compile_uf(force_type t);

/**
* @brief compiles the body update program
*/
template <typename T>
func_ptr<T> compile_ub(update_type t);

template <typename T> 
void upload_force_kernel(cl_kernel& k, Simulation<T>* _s);

template <typename T> 
void upload_update_kernel(cl_kernel& k, Simulation<T>* _s, double c = 1, double d = 1);

/**
* @brief: gets_the maximum amount of data we can put on the GPU (assuming everything is a uint32_t)
* @param device: device to scan
* @returns the amount of uint32_t to load
*/
static size_t compute_max_chunk_size(cl_device_id device);

/**
* @brief merges two already sorted arrays
* @param left first array
* @param L first array size
* @param right second array
* @param R second array size
* @param merged ptr to the merged array
*/
static void merge_arrays(const uint32_t* left, size_t L,const uint32_t* right, size_t R,uint32_t* merged);

cl_program compile_sorting_kernels();

/**
* @brief sorts an array using the gpu
* @param input: ptr to the start of the arrya
* @param output: array onto which to write
* @param N: size of the array
*/
void sort(uint64_t* input, uint64_t* output, size_t N);
}


}

#include "Aster/impl/tc_impl_cl.tpp"
#include "Aster/impl/tool_chain_impl.tpp"
#include "Aster/impl/SABA.tpp"