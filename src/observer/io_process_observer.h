#ifndef UG_PLUGIN_XBRAIDFORUG4_IO_PROCESS_OBSERVER_H
#define UG_PLUGIN_XBRAIDFORUG4_IO_PROCESS_OBSERVER_H

#include <map>
#include <tuple>

#include "lib_disc/function_spaces/grid_function.h"

#include "../interface/observer_xbraid.h"
#include "../util/parallel_io_gridfunction.h"


namespace ug { namespace XBraidForUG4 {

    template<typename TDomain, typename TAlgebra>
     class IO_ProcessObserver : public IXBraidTimeIntegratorObserver<TDomain, TAlgebra> {
    public:

        //--------------------------------------------------------------------------------------------------------------

        typedef GridFunction<TDomain, TAlgebra> T_GridFunction;
        typedef SmartPtr<T_GridFunction> SP_GridFunction;

        typedef PIOGridFunction<TDomain, TAlgebra> T_IOOutput;

        typedef std::tuple<int, int, int> TKey;

        //--------------------------------------------------------------------------------------------------------------

        const char *m_filename;

        std::map<TKey, int> map;

        //--------------------------------------------------------------------------------------------------------------

        explicit IO_ProcessObserver(const char *filename) : IXBraidTimeIntegratorObserver<TDomain, TAlgebra>() {
            m_filename = filename;
        }
        ~IO_ProcessObserver() override = default;

        //--------------------------------------------------------------------------------------------------------------

        void set_filename(const char *filename) {
            this->m_filename = filename;
        }

        bool step_process(SP_GridFunction u, int index, number time, number dt) override {
            T_IOOutput io = T_IOOutput();
            std::stringstream ss;
            ss << m_filename << "_t" << index << ".gridfunction";
            io.write(u, ss.str().c_str());
            return true;
        };

        bool step_process(SP_GridFunction u, int index, double time, double dt, int iteration, int level) override {

            int count = 0; // todo check if this is still neccessary for the newer xbraid version
            auto tuple = std::make_tuple(index, iteration, level);
            auto it = map.find(tuple);
            if (it != map.end()) {
                count = it->second;
                count += 1;
                map[tuple] = count;
            } else {
                count = 0;
                map.emplace(tuple, 0);
            }

            std::stringstream ss;
            ss << m_filename << "_k" << iteration << "_l" << level << "_c" << count << "_t" << index
               << ".gridfunction";
            T_IOOutput io = T_IOOutput();
            io.write(u, ss.str().c_str());

            return true;
        };
        //--------------------------------------------------------------------------------------------------------------
    };
}}
#endif //UG_PLUGIN_XBRAIDFORUG4_IO_PROCESS_OBSERVER_H
