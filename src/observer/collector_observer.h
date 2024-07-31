#ifndef UG_PLUGIN_XBRAIDFORUG4_COLLECTOR_OBSERVER_H
#define UG_PLUGIN_XBRAIDFORUG4_COLLECTOR_OBSERVER_H

#include <vector>

#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/time_disc/time_integrator_observers/time_integrator_observer_interface.h"



namespace ug { namespace XBraidForUG4 {

    template <typename TDomain, typename TAlgebra>
    class TimeIntegratorObserverCollector : public ITimeIntegratorObserver<TDomain, TAlgebra> {
    public:

        //--------------------------------------------------------------------------------------------------------------

        typedef GridFunction<TDomain, TAlgebra> T_GridFunction;
        typedef SmartPtr<T_GridFunction> SP_GridFunction;


        typedef ITimeIntegratorObserver<TDomain, TAlgebra> T_ITimeIntegratorObserver;
        typedef SmartPtr<T_ITimeIntegratorObserver> SP_ITimeIntegratorObserver;

        //--------------------------------------------------------------------------------------------------------------

        std::vector<SP_ITimeIntegratorObserver> lst;

        //--------------------------------------------------------------------------------------------------------------

        TimeIntegratorObserverCollector() = default;
        ~TimeIntegratorObserverCollector() override = default;

        //--------------------------------------------------------------------------------------------------------------

        void attach_observer(SP_ITimeIntegratorObserver observer) {
            this->lst.emplace_back(observer);
        }


        bool step_process(SP_GridFunction u, int index, double time, double dt) override {
            for (size_t k = 0; k < this->lst.size(); k++) {
                this->lst[k]->step_process(u, index, time, dt);
            }
            return true;
        };

        //--------------------------------------------------------------------------------------------------------------
    };
}}
#endif //UG_PLUGIN_XBRAIDFORUG4_COLLECTOR_OBSERVER_H
