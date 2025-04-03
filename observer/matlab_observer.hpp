#ifndef UGPLUGIN_XBRAIDFORUG4_OBSERVER_MATLAB_OBSERVER_HPP
#define UGPLUGIN_XBRAIDFORUG4_OBSERVER_MATLAB_OBSERVER_HPP


#include "lib_disc/time_disc/time_integrator_observers/time_integrator_observer_interface.h"




namespace ug{ namespace xbraid {

    template<typename TDomain, typename TAlgebra>
    class MATLAB_Observer : public ITimeIntegratorObserver<TDomain, TAlgebra> {
    public:

        //--------------------------------------------------------------------------------------------------------------

        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;


        using T_ParallelLogger = ParallelLogger ;
        using SP_ParallelLogger = SmartPtr<T_ParallelLogger> ;


        //--------------------------------------------------------------------------------------------------------------

        explicit MATLAB_Observer(SP_ParallelLogger log) : ITimeIntegratorObserver<TDomain, TAlgebra>(){
            this->logger_ = log;
        }

        ~MATLAB_Observer() override = default;

        //--------------------------------------------------------------------------------------------------------------



        bool step_process(SP_GridFunction u, int index, double time, double dt) override {
            int sz = u->size();
            this->logger_->o << "u_" << index << " = [";
            for (int i = 0; i < sz - 1; i++) {
                this->logger_->o << std::setw(7) << (*u)[i] << ", ";
            }
            this->logger_->o << std::setw(7) << (*u)[sz - 1];
            this->logger_->o << "]" << std::endl;
            return true;
        };

        //--------------------------------------------------------------------------------------------------------------

        SP_ParallelLogger logger_;

        //--------------------------------------------------------------------------------------------------------------
    };
}}
#endif