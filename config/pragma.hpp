#ifndef UGPLUGIN_XBRAIDFORUG4_CONFIG_PRAGMA_HPP
#define UGPLUGIN_XBRAIDFORUG4_CONFIG_PRAGMA_HPP

#include "compile_settings.hpp"
#include "execution_settings.hpp"

#ifdef MDEBUG
    #define __debug(expression) {expression}
#else
    #define __debug(expression)
#endif



#include <iostream>

// #define MDEBUG
//#define WRITE_SCRIPT
//#define TRACE_RECVTIME 1

// #define ExtendedInfo
// #define SEND_RECV_TIMES
#define FEATURE_SPATIAL_REFINE

#ifdef ExtendedInfo
    #define extended_info(expression) expression
#else
    #define extended_info(expression)
#endif


#ifdef WRITE_SCRIPT
    #define write_script(expression) expression
#else
    #define write_script(expression)
#endif


#ifdef SEND_RECV_TIMES
#define __send_recv_times(expression) expression
#else
#define __send_recv_times(expression)
#endif

#ifdef FEATURE_SPATIAL_REFINE
#else
#endif



#define __debug_print(expression) {std::cout << expression << std::endl << std::flush}





#endif