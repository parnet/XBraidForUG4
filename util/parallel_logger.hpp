#ifndef UGPLUGIN_XBRAIDFORUG4_UTIL_PARALLEL_LOGGER_HPP
#define UGPLUGIN_XBRAIDFORUG4_UTIL_PARALLEL_LOGGER_HPP

#include <fstream>

#include "common/util/smart_pointer.h"

#include "core/space_time_communicator.hpp"

namespace ug { namespace xbraid {

    class ParallelLogger {
    public:
        //--------------------------------------------------------------------------------------------------------------

        using SP_SpaceTimeCommunicator = SmartPtr<SpaceTimeCommunicator> ;




        //--------------------------------------------------------------------------------------------------------------

        ParallelLogger() = default;
        ~ParallelLogger() = default;

        //--------------------------------------------------------------------------------------------------------------

        void set_comm(SP_SpaceTimeCommunicator comm) {
            comm_set_ = true;
            this->comm_ = comm;
        }

        void set_filename(const char *filename) {
            this->filename_ = filename;
        }

        void init() {
            if (this->init_ == 0) {
                std::stringstream ss;
                if(comm_set_) {
                    ss << this->filename_ << "_" << this->comm_->get_temporal_rank() << ".output";
                } else {
                    ss << this->filename_ << "_" << 0 << ".output";
                }
                this->o.open(ss.str());
            }
            this->init_++;
        }

        void release() {
            this->init_--;
            if (this->init_ == 0) {
                this->o << std::endl << std::endl
                        << "finished" << std::endl << std::flush;
            }
        }

        void write(const char *content) {
            this->o << content << std::endl;
        }
        //--------------------------------------------------------------------------------------------------------------

        SP_SpaceTimeCommunicator comm_;

        bool comm_set_ = false;
        int init_ = 0;
        const char *filename_ = "job";

        std::ofstream o;
    };
}}
#endif