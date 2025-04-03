#ifndef UGPLUGIN_XBRAIDFORUG4_CORE_DISCARD_STANDARD_STREAM_HPP
#define UGPLUGIN_XBRAIDFORUG4_CORE_DISCARD_STANDARD_STREAM_HPP

#include "space_time_communicator.hpp"





namespace ug{ namespace xbraid {


    class DiscardStandardStream {
    public:
        //--------------------------------------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------------------------------------


        std::streambuf* stdcout_;

        bool active_ = false;

        //--------------------------------------------------------------------------------------------------------------

        DiscardStandardStream() = default;

        ~DiscardStandardStream() = default;

        //--------------------------------------------------------------------------------------------------------------


        void apply() {
            if (active_) { return; }

            this->stdcout_ = std::cout.rdbuf(); // save state
            std::cout.rdbuf(nullptr); // replace state

            this->active_ = true;

        }

        void undo() {
            if(!this->active_){return;}

            std::cout.rdbuf(this->stdcout_); // replace state

            this->active_ = false;
        }

        //--------------------------------------------------------------------------------------------------------------
    };

}}
#endif