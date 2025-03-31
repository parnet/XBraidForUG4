#ifndef UGPLUGIN_XBRAIDFORUG4_OBSERVER_IO_OBSERVER_HPP
#define UGPLUGIN_XBRAIDFORUG4_OBSERVER_IO_OBSERVER_HPP

#include <tuple>
#include "lib_disc/time_disc/time_integrator_observers/time_integrator_observer_interface.h"
#include "util/parallel_io_gridfunction.hpp"



namespace ug{ namespace xbraid {

    template<typename TDomain, typename TAlgebra>
    class IO_Observer : public ITimeIntegratorObserver<TDomain, TAlgebra> {
    public:

        //--------------------------------------------------------------------------------------------------------------

        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;

        using T_IOOutput = PIOGridFunction<TDomain, TAlgebra> ;

        using T_Key = std::tuple<int, int, int> ;

        //--------------------------------------------------------------------------------------------------------------

        const char *m_filename;

        //--------------------------------------------------------------------------------------------------------------

        explicit IO_Observer(const char *filename) : ITimeIntegratorObserver<TDomain, TAlgebra>() {
            m_filename = filename;
        }


        ~IO_Observer() override = default;

        //--------------------------------------------------------------------------------------------------------------

        void set_filename(const char *filename) {
            this->m_filename = filename;
        }

        bool step_process(SP_GridFunction u, int index, double time, double dt) override {
            T_IOOutput io = T_IOOutput();
            std::stringstream ss;
            ss << m_filename << "_t" << index << ".gridfunction";
            io.write(u, ss.str().c_str());
            return true;
        }

        //--------------------------------------------------------------------------------------------------------------
    };
}}
#endif