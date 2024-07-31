#ifndef UG_PLUGIN_XBRAIDFORUG4_IO_OBSERVER_H
#define UG_PLUGIN_XBRAIDFORUG4_IO_OBSERVER_H


#include <tuple>
#include "lib_disc/time_disc/time_integrator_observers/time_integrator_observer_interface.h"
#include "../util/parallel_io_gridfunction.h"


namespace ug {namespace XBraidForUG4 {

    template<typename TDomain, typename TAlgebra>
    class IO_Observer : public ITimeIntegratorObserver<TDomain, TAlgebra> {
    public:

        //--------------------------------------------------------------------------------------------------------------

        typedef GridFunction<TDomain, TAlgebra> T_GridFunction;
        typedef SmartPtr<T_GridFunction> SP_GridFunction;

        typedef PIOGridFunction<TDomain, TAlgebra> T_IOOutput;

        typedef std::tuple<int, int, int> T_Key;

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
#endif //UG_PLUGIN_XBRAIDFORUG4_IO_OBSERVER_H
