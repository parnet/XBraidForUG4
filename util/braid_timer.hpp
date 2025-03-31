#ifndef UGPLUGIN_XBRAIDFORUG4_UTIL_BRAID_TIMER_HPP
#define UGPLUGIN_XBRAIDFORUG4_UTIL_BRAID_TIMER_HPP

#include <chrono>





namespace ug { namespace xbraid {

    class BraidTimer {
    public:
        //--------------------------------------------------------------------------------------------------------------

        std::chrono::high_resolution_clock::time_point t0;
        std::chrono::high_resolution_clock::time_point t1;
        std::chrono::high_resolution_clock::time_point t2;

        //--------------------------------------------------------------------------------------------------------------

        BraidTimer() {
            t0 = std::chrono::high_resolution_clock::now();
            t1 = std::chrono::high_resolution_clock::now();
            t2 = std::chrono::high_resolution_clock::now();
        }

        ~BraidTimer() = default;

        //--------------------------------------------------------------------------------------------------------------

        void start() {
            t0 = std::chrono::high_resolution_clock::now();
        }

        void stop() {
            t1 = std::chrono::high_resolution_clock::now();
        }

        void now(double& total, double& diff) {
            t2 = std::chrono::high_resolution_clock::now();
            const std::chrono::duration<double> t0diff = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t0);
            total = t0diff.count();
            diff = 0;
        }


        double get() const {
            const std::chrono::duration<double> difference = std::chrono::duration_cast<std::chrono::duration<double> >(t1 - t0);
            return difference.count();
        }

        //--------------------------------------------------------------------------------------------------------------
    };

}}

#endif