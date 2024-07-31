#ifndef REPLACE_STANDARD_STREAM_H
#define REPLACE_STANDARD_STREAM_H

#include "space_time_communicator.h"

namespace ug { namespace XBraidForUG4 {


class ReplaceStandardStream {
public:
    //--------------------------------------------------------------------------------------------------------------

    typedef SmartPtr<SpaceTimeCommunicator> SP_SpaceTimeCommunicator;

    //--------------------------------------------------------------------------------------------------------------

    SP_SpaceTimeCommunicator m_comm;

    std::streambuf* stdcout;
    std::ofstream subbuff; // todo close file -> undo or destructor

    //--------------------------------------------------------------------------------------------------------------

    ReplaceStandardStream() = default;

    ~ReplaceStandardStream() = default;

    //--------------------------------------------------------------------------------------------------------------

    void set_space_time_comm(SP_SpaceTimeCommunicator comm) {
        this->m_comm = comm;
    }

    void apply() {
        std::stringstream ss;
        ss << "std_output_" << this->m_comm->get_temporal_rank() << ".cout";
        this->subbuff.open(ss.str());

        this->stdcout = std::cout.rdbuf(); // save state
        std::cout.rdbuf(subbuff.rdbuf()); // replace state

    }

    void undo() const {
        std::cout.rdbuf(this->stdcout); // replace state
    }

    //--------------------------------------------------------------------------------------------------------------
};
}}
#endif //REPLACE_STANDARD_STREAM_H
