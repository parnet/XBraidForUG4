#ifndef UGPLUGIN_XBRAIDFORUG4_CORE_REPLACE_STANDARD_STREAM_H
#define UGPLUGIN_XBRAIDFORUG4_CORE_REPLACE_STANDARD_STREAM_H

#include "space_time_communicator.hpp"
//2025-03 #include "discard_standard_stream.h"




namespace ug{ namespace xbraid {

class ReplaceStandardStream {
public:
    //--------------------------------------------------------------------------------------------------------------

    using SP_SpaceTimeCommunicator = SmartPtr<SpaceTimeCommunicator> ;

    //--------------------------------------------------------------------------------------------------------------

    SP_SpaceTimeCommunicator comm_;

    std::streambuf* stdcout_;
    std::ofstream subbuff_;

    //--------------------------------------------------------------------------------------------------------------

    ReplaceStandardStream() = default;

    ~ReplaceStandardStream() = default;

    //--------------------------------------------------------------------------------------------------------------

    void set_space_time_comm(SP_SpaceTimeCommunicator comm) {
        this->comm_ = comm;
    }

    void apply() {
        std::stringstream ss;
        ss << "std_output_" << this->comm_->get_temporal_rank() << ".cout";
        this->subbuff_.open(ss.str());

        this->stdcout_ = std::cout.rdbuf(); // save state
        std::cout.rdbuf(subbuff_.rdbuf()); // replace state

    }

    void undo() const {
        std::cout.rdbuf(this->stdcout_); // replace state
    }

    //--------------------------------------------------------------------------------------------------------------
};

}}
#endif
