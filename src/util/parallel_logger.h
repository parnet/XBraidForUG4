#ifndef UG_PLUGIN_XBRAIDFORUG4_PARALOG_H
#define UG_PLUGIN_XBRAIDFORUG4_PARALOG_H

#include <iostream>
#include <fstream>

#include "common/util/smart_pointer.h"

#include "../core/space_time_communicator.h"




namespace ug { namespace XBraidForUG4 {

    class ParallelLogger {
    public:
        //--------------------------------------------------------------------------------------------------------------

        typedef SmartPtr<SpaceTimeCommunicator> SP_SpaceTimeCommunicator;

        //--------------------------------------------------------------------------------------------------------------

        SP_SpaceTimeCommunicator m_comm;

        bool m_comm_set = false;
        int m_init = 0;
        const char *m_filename = "job";

        std::ofstream o;


        //--------------------------------------------------------------------------------------------------------------

        ParallelLogger() = default;
        ~ParallelLogger() = default;

        //--------------------------------------------------------------------------------------------------------------

        void set_comm(SP_SpaceTimeCommunicator comm) {
            m_comm_set = true;
            this->m_comm = comm;
        }

        void set_filename(const char *filename) {
            this->m_filename = filename;
        }

        void init() {
            if (this->m_init == 0) {
                std::stringstream ss;
                if(m_comm_set) {
                    ss << this->m_filename << "_" << this->m_comm->get_temporal_rank() << ".output";
                } else {
                    ss << this->m_filename << "_" << 0 << ".output";
                }
                this->o.open(ss.str());
            }
            this->m_init++;
        }

        void release() {
            this->m_init--;
            if (this->m_init == 0) {
                this->o << std::endl << std::endl << "finished" << std::endl << std::flush;
            }
        }

        void write(const char *content) {
            this->o << content << std::endl;
        }

    };
}}
#endif //UG_PLUGIN_XBRAIDFORUG4_PARALOG_H
