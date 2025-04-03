#ifndef UGPLUGIN_XBRAIDFORUG4_POROELASTICITY_BRAID_BIOT_CONTROL_HPP
#define UGPLUGIN_XBRAIDFORUG4_POROELASTICITY_BRAID_BIOT_CONTROL_HPP

#include "braid_biot_precomputed.hpp"
#include "Poroelasticity/biot_tools.h"
#include "Poroelasticity/barry_mercer.h"

#include "interface/observer_xbraid.hpp"

namespace ug { namespace xbraid { namespace poro {

    /**
     * A specialized observer class to write out the result of the BarryMercerProblem as a parallel io gridfunction
     * @tparam TDomain
     * @tparam TAlgebra
     */
    template <typename TDomain, typename TAlgebra>
    class BraidBiotCheck final : public IXBraidTimeIntegratorObserver<TDomain, TAlgebra> {
    public:

        //--------------------------------------------------------------------------------------------------------------

        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;

        using T_Problem = Poroelasticity::BarryMercerProblem<TDomain, TAlgebra> ;
        using SP_Problem = SmartPtr<T_Problem> ;

        using T_PIOGridFunction = PIOGridFunction<TDomain, TAlgebra> ;
        using SP_PIOGridFunction = SmartPtr<T_PIOGridFunction> ;


        //--------------------------------------------------------------------------------------------------------------

        BraidBiotCheck() : IXBraidTimeIntegratorObserver<TDomain, TAlgebra>() {};

        ~BraidBiotCheck() override = default;

        //--------------------------------------------------------------------------------------------------------------

        void set_problem(SP_Problem problem) {
            this->problem_ = problem;
        }

        /*void set_napprox(int napprox) {
            this->napprox = napprox;
        }*/

        void set_filename(const char * filename) {
            this->filename_ = filename;
        }

        virtual bool step_process(SP_GridFunction u, int index, double time, double dt) {
            //m_problem->m_errData.napprox = this->napprox;

            //æ SP_GridFunction solution = m_problem->compute_solution(u, index, time);
            //                 m_problem->m_errData.napprox = this->napprox;
            //            m_problem->m_errData.iteration = -1;
            //            m_problem->m_errData.level = -1;
            //            m_problem->post_processing(u, index, time);
            SP_GridFunction solution; // = m_problem->compute_solution(u, index, time);

            std::stringstream ss;
            ss << this->filename_<<"_"<<index;
            SP_PIOGridFunction pio = SP_PIOGridFunction();
            pio->write(solution,ss.str().c_str());

            return false;
        };

        bool step_process(SP_GridFunction u, int index, double time, double dt, int interation, int level) override {
            this->step_process(u,index,time,dt);
            //            m_problem->m_errData.napprox = this->napprox;
            //            m_problem->m_errData.iteration = iteration;
            //            m_problem->m_errData.level = level;
            //
            //            m_problem->post_processing(u, index, time);
            //
            //            m_problem->m_errData.iteration = -1;
            //            m_problem->m_errData.level = -1;
            return false;
        };

        SP_GridFunction lua_write(SP_GridFunction u, int index, double time, double dt) {
            //m_problem->m_errData.napprox = this->napprox;

            //æ SP_GridFunction solution = m_problem->compute_solution(u, index, time);
            SP_GridFunction solution; // = m_problem->compute_solution(u, index, time);

            std::stringstream ss;
            ss << this->filename_<<"_"<<index;
            SP_PIOGridFunction pio = SP_PIOGridFunction();
            pio->write(solution,ss.str().c_str());

            return solution;
        }

        //--------------------------------------------------------------------------------------------------------------

        SP_Problem problem_;
        int napprox_ = 16;
        const char * filename_ = "BarryMercer_2D";
        //--------------------------------------------------------------------------------------------------------------
    };
}}}


#endif