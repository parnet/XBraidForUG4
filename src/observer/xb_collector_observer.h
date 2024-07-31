#ifndef UG_PLUGIN_XBRAIDFORUG4_XB_COLLECTOR_PROCESS_OBSERVER_H
#define UG_PLUGIN_XBRAIDFORUG4_XB_COLLECTOR_PROCESS_OBSERVER_H

#include "../interface/observer_xbraid.h"


namespace ug { namespace XBraidForUG4 {
    template <typename TDomain, typename TAlgebra>
    class XBraidTimeIntegratorObserverCollector : public IXBraidTimeIntegratorObserver<TDomain, TAlgebra> {
    public:

        //--------------------------------------------------------------------------------------------------------------

        typedef GridFunction<TDomain, TAlgebra> T_GridFunction;
        typedef SmartPtr<T_GridFunction> SP_GridFunction;

        typedef ITimeIntegratorObserver<TDomain, TAlgebra> T_TimeIntegratorObserver;
        typedef SmartPtr<T_TimeIntegratorObserver> SP_TimeIntegratorObserver;

        typedef IXBraidTimeIntegratorObserver<TDomain, TAlgebra> T_XBraidTimeIntegratorObserver;
        typedef SmartPtr<T_XBraidTimeIntegratorObserver> SP_XBraidTimeIntegratorObserver;

        //--------------------------------------------------------------------------------------------------------------

        std::vector<SP_TimeIntegratorObserver> lst;
        std::vector<SP_XBraidTimeIntegratorObserver> xb_lst;

        //--------------------------------------------------------------------------------------------------------------

        XBraidTimeIntegratorObserverCollector() {};

        ~XBraidTimeIntegratorObserverCollector() override {};

        //--------------------------------------------------------------------------------------------------------------

        void attach_observer(SP_XBraidTimeIntegratorObserver observer) {
            this->xb_lst.emplace_back(observer);
        }

        void attach_common_observer(SP_TimeIntegratorObserver observer) {
            this->lst.emplace_back(observer);
        }


        bool step_process(SP_GridFunction u, int index, number time, number dt) override {
            for (size_t k = 0; k < this->lst.size(); k++) {
                this->lst[k]->step_process(u, index, time,dt);
            }

            for (size_t k = 0; k < this->xb_lst.size(); k++) {
                this->xb_lst[k]->step_process(u, index, time,dt);
            }

            return true;
        };

        bool step_process(SP_GridFunction u, int index, number time, number dt , int iteration, int level) override {
            for (size_t k = 0; k < this->xb_lst.size(); k++) {
                this->xb_lst[k]->step_process(u, index, time, dt, iteration, level);
            }
            return true;
        };

        //--------------------------------------------------------------------------------------------------------------
    };
}}
#endif //UG_PLUGIN_XBRAIDFORUG4_XB_COLLECTOR_PROCESS_OBSERVER_H
