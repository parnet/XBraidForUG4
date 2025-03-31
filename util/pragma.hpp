#ifndef UGPLUGIN_XBRAIDFORUG4_COMMON_PRAGMA_HPP
#define UGPLUGIN_XBRAIDFORUG4_COMMON_PRAGMA_HPP


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


#ifdef MDEBUG
	#define __debug(expression) {expression}
#else
	#define __debug(expression)
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





#define TRACE_TIMINGS 1
#if TRACE_TIMINGS == 1
	#define StartOperationTimer(opt) this->m_time_log_manager.get(opt).start()
	#define StopOperationTimer(opt) this->m_time_log_manager.get(opt).stop()
	#define StartLevelOperationTimer(opt, l) this->m_time_log_manager.get(opt,l).start()
	#define StopLevelOperationTimer(opt, l) this->m_time_log_manager.get(opt,l).stop()
#else
	#define StartOperationTimer(opt)
	#define StopOperationTimer(opt)
	#define StartLevelOperationTimer(opt, l)
	#define StopLevelOperationTimer(opt, l)
#endif




#endif