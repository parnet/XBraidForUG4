#ifndef UGPLUGIN_XBRAIDFORUG4_CONFIG_EXECUTION_SETTINGS_HPP
#define UGPLUGIN_XBRAIDFORUG4_CONFIG_EXECUTION_SETTINGS_HPP




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