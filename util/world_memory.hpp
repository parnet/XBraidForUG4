#ifndef UGPLUGIN_XBRAIDFORUG4_UTIL_WORLD_MEMORY_H
#define UGPLUGIN_XBRAIDFORUG4_UTIL_WORLD_MEMORY_H







namespace ug { namespace xbraid {

        unsigned long get_world_memory_consumed();
        unsigned long get_world_memory_peak();
        unsigned long get_spatial_memory_consumed();
        void get_world_memory_distribution();
        void get_spatial_memory_distribution();

}}

#endif