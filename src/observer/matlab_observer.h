#ifndef UG_PLUGIN_XBRAIDFORUG4_MATLAB_OBSERVER_H
#define UG_PLUGIN_XBRAIDFORUG4_MATLAB_OBSERVER_H


#include "lib_disc/time_disc/time_integrator_observers/time_integrator_observer_interface.h"

namespace ug { namespace XBraidForUG4 {
    template<typename TDomain, typename TAlgebra>
    class MATLAB_Observer : public ITimeIntegratorObserver<TDomain, TAlgebra> {
    public:

        //--------------------------------------------------------------------------------------------------------------

        typedef GridFunction<TDomain, TAlgebra> T_GridFunction;
        typedef SmartPtr<T_GridFunction> SP_GridFunction;


        typedef ParallelLogger T_ParallelLogger;
        typedef SmartPtr<T_ParallelLogger> SP_ParallelLogger;

        //--------------------------------------------------------------------------------------------------------------

        SP_ParallelLogger logger;

        //--------------------------------------------------------------------------------------------------------------

        explicit MATLAB_Observer(SP_ParallelLogger log) : ITimeIntegratorObserver<TDomain, TAlgebra>(){
            this->logger = log;
        }

        ~MATLAB_Observer() override = default;

        //--------------------------------------------------------------------------------------------------------------



        bool step_process(SP_GridFunction u, int index, double time, double dt) override {
            int sz = u->size();
            this->logger->o << "u_" << index << " = [";
            for (int i = 0; i < sz - 1; i++) {
                this->logger->o << std::setw(7) << (*u)[i] << ", ";
            }
            this->logger->o << std::setw(7) << (*u)[sz - 1];
            this->logger->o << "]" << std::endl;
            return true;
        };

        //--------------------------------------------------------------------------------------------------------------
    };
}}
#endif //UG_PLUGIN_XBRAIDFORUG4_MATLAB_OBSERVER_H
