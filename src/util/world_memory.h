#ifndef UG_PLUGIN_XBRAIDFORUG4_WORLD_MEMORY_H
#define UG_PLUGIN_XBRAIDFORUG4_WORLD_MEMORY_H



namespace ug { namespace XBraidForUG4 {

        unsigned long get_world_memory_consumed();
        unsigned long get_world_memory_peak();
        unsigned long get_spatial_memory_consumed();
        void get_world_memory_distribution();
        void get_spatial_memory_distribution();

}}

#endif //UG_PLUGIN_XBRAIDFORUG4_WORLD_MEMORY_H
