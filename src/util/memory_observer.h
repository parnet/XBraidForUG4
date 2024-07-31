#ifndef UG_PLUGIN_XBRAIDFORUG4_MEMORY_OBSERVER_H
#define UG_PLUGIN_XBRAIDFORUG4_MEMORY_OBSERVER_H



namespace ug { namespace XBraidForUG4 {

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

#endif //UG_PLUGIN_XBRAIDFORUG4_MEMORY_OBSERVER_H
