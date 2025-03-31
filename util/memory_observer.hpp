#ifndef UGPLUGIN_XBRAIDFORUG4_UTIL_MEMORY_OBSERVER_HPP
#define UGPLUGIN_XBRAIDFORUG4_UTIL_MEMORY_OBSERVER_HPP







namespace ug { namespace xbraid {

    unsigned long parseLine(char* line);
    unsigned long get_virtual_memory_total();
    unsigned long get_virtual_memory_used();
    unsigned long get_virtual_memory_consumed();
    unsigned long get_virtual_memory_peak();
    unsigned long get_physical_memory_total();
    unsigned long get_physical_memory_used();
    unsigned long get_physical_memory_consumed();
    unsigned long get_physical_memory_peak();

}}

#endif