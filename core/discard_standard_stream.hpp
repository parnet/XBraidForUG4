#ifndef UGPLUGIN_XBRAIDFORUG4_CORE_DISCARD_STANDARD_STREAM_HPP
#define UGPLUGIN_XBRAIDFORUG4_CORE_DISCARD_STANDARD_STREAM_HPP

#include "space_time_communicator.hpp"





namespace ug{ namespace xbraid {


    class DiscardStandardStream {
    public:
        //--------------------------------------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------------------------------------


        std::streambuf* stdcout;

        bool active = false;

        //--------------------------------------------------------------------------------------------------------------

        DiscardStandardStream() = default;

        ~DiscardStandardStream() = default;

        //--------------------------------------------------------------------------------------------------------------


        void apply() {
            if (this->active) { return; }

            this->stdcout = std::cout.rdbuf(); // save state
            std::cout.rdbuf(nullptr); // replace state

            this->active = true;

        }

        void undo() {
            if(!this->active){return;}

            std::cout.rdbuf(this->stdcout); // replace state

            this->active = false;
        }

        //--------------------------------------------------------------------------------------------------------------
    };

}}
#endif