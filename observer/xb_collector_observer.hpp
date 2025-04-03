#ifndef UGPLUGIN_XBRAIDFORUG4_OBSERVER_XB_COLLECTOR_PROCESS_OBSERVER_HPP
#define UGPLUGIN_XBRAIDFORUG4_OBSERVER_XB_COLLECTOR_PROCESS_OBSERVER_HPP

#include "interface/observer_xbraid.hpp"


namespace ug{ namespace xbraid {
    template <typename TDomain, typename TAlgebra>
    class XBraidTimeIntegratorObserverCollector : public IXBraidTimeIntegratorObserver<TDomain, TAlgebra> {
    public:

        //--------------------------------------------------------------------------------------------------------------

        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;

        using T_TimeIntegratorObserver = ITimeIntegratorObserver<TDomain, TAlgebra> ;
        using SP_TimeIntegratorObserver= SmartPtr<T_TimeIntegratorObserver> ;

        using T_XBraidTimeIntegratorObserver = IXBraidTimeIntegratorObserver<TDomain, TAlgebra> ;
        using SP_XBraidTimeIntegratorObserver = SmartPtr<T_XBraidTimeIntegratorObserver> ;

        //--------------------------------------------------------------------------------------------------------------

        XBraidTimeIntegratorObserverCollector() = default;

        ~XBraidTimeIntegratorObserverCollector() override = default;

        //--------------------------------------------------------------------------------------------------------------

        void attach_observer(SP_XBraidTimeIntegratorObserver observer) {
            this->xb_lst_.emplace_back(observer);
        }

        void attach_common_observer(SP_TimeIntegratorObserver observer) {
            this->lst_.emplace_back(observer);
        }


        bool step_process(SP_GridFunction u, int index, number time, number dt) override {
            for (size_t k = 0; k < this->lst_.size(); k++) {
                this->lst_[k]->step_process(u, index, time,dt);
            }

            for (size_t k = 0; k < this->xb_lst_.size(); k++) {
                this->xb_lst_[k]->step_process(u, index, time,dt);
            }

            return true;
        };

        bool step_process(SP_GridFunction u, int index, number time, number dt , int iteration, int level) override {
            for (size_t k = 0; k < this->xb_lst_.size(); k++) {
                this->xb_lst_[k]->step_process(u, index, time, dt, iteration, level);
            }
            return true;
        };

        //--------------------------------------------------------------------------------------------------------------

        std::vector<SP_TimeIntegratorObserver> lst_;
        std::vector<SP_XBraidTimeIntegratorObserver> xb_lst_;

        //--------------------------------------------------------------------------------------------------------------
    };
}}
#endif