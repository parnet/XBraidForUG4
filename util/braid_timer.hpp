#ifndef UGPLUGIN_XBRAIDFORUG4_UTIL_BRAID_TIMER_HPP
#define UGPLUGIN_XBRAIDFORUG4_UTIL_BRAID_TIMER_HPP

#include <chrono>





namespace ug { namespace xbraid {

    class BraidTimer {
    public:

        //--------------------------------------------------------------------------------------------------------------

        BraidTimer() {
            t0_ = std::chrono::high_resolution_clock::now();
            t1_ = std::chrono::high_resolution_clock::now();
            t2_ = std::chrono::high_resolution_clock::now();
        }

        ~BraidTimer() = default;

        //--------------------------------------------------------------------------------------------------------------

        void start() {
            t0_ = std::chrono::high_resolution_clock::now();
        }

        void stop() {
            t1_ = std::chrono::high_resolution_clock::now();
        }

        void now(double& total, double& diff) {
            t2_ = std::chrono::high_resolution_clock::now();
            const std::chrono::duration<double> t0diff = std::chrono::duration_cast<std::chrono::duration<double>>(t2_ - t0_);
            total = t0diff.count();
            diff = 0; // todo set value for difference from last time to now or remove second parameter and create return value
        }


        double get() const {
            const std::chrono::duration<double> difference = std::chrono::duration_cast<std::chrono::duration<double> >(t1_ - t0_);
            return difference.count();
        }
        //--------------------------------------------------------------------------------------------------------------

        std::chrono::high_resolution_clock::time_point t0_;
        std::chrono::high_resolution_clock::time_point t1_;
        std::chrono::high_resolution_clock::time_point t2_;

        //--------------------------------------------------------------------------------------------------------------
    };

}}

#endif