#ifndef UGPLUGIN_XBRAIDFORUG4_IO_PROCESS_OBSERVER_HPP
#define UGPLUGIN_XBRAIDFORUG4_IO_PROCESS_OBSERVER_HPP

#include <map>
#include <tuple>

#include "lib_disc/function_spaces/grid_function.h"

#include "interface/observer_xbraid.hpp"
#include "util/parallel_io_gridfunction.hpp"

namespace ug{ namespace xbraid {

    template<typename TDomain, typename TAlgebra>
     class IO_ProcessObserver : public IXBraidTimeIntegratorObserver<TDomain, TAlgebra> {
    public:

        //--------------------------------------------------------------------------------------------------------------

        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;

        using T_IOOutput = PIOGridFunction<TDomain, TAlgebra> ;

        using TKey = std::tuple<int, int, int> ;



        //--------------------------------------------------------------------------------------------------------------

        explicit IO_ProcessObserver(const char *filename) : IXBraidTimeIntegratorObserver<TDomain, TAlgebra>() {
            filename_ = filename;
        }
        ~IO_ProcessObserver() override = default;

        //--------------------------------------------------------------------------------------------------------------

        void set_filename(const char *filename) {
            this->filename_ = filename;
        }

        bool step_process(SP_GridFunction u, int index, number time, number dt) override {
            T_IOOutput io = T_IOOutput();
            std::stringstream ss;
            ss << filename_ << "_t" << index << ".gridfunction";
            io.write(u, ss.str().c_str());
            return true;
        };

        bool step_process(SP_GridFunction u, int index, double time, double dt, int iteration, int level) override {

            int count = 0;
            auto tuple = std::make_tuple(index, iteration, level);
            auto it = map_.find(tuple);
            if (it != map_.end()) {
                count = it->second;
                count += 1;
                map_[tuple] = count;
            } else {
                count = 0;
                map_.emplace(tuple, 0);
            }

            std::stringstream ss;
            ss << filename_ << "_k" << iteration << "_l" << level << "_c" << count << "_t" << index
               << ".gridfunction";
            T_IOOutput io = T_IOOutput();
            io.write(u, ss.str().c_str());

            return true;
        };
        //--------------------------------------------------------------------------------------------------------------

        const char *filename_;

        std::map<TKey, int> map_;
        //--------------------------------------------------------------------------------------------------------------
    };
}}
#endif