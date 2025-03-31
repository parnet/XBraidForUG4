#ifndef UGPLUGIN_XBRAIDFORUG4_DRIVER_BASIC_DRIVER_HPP
#define UGPLUGIN_XBRAIDFORUG4_DRIVER_BASIC_DRIVER_HPP



#include "interface/residual_timeintegrator.hpp"
#include "gridfunction_base.hpp"



namespace ug{ namespace xbraid {

    template <typename TDomain, typename TAlgebra>
    class BasicDriver : public BraidGridFunctionBase<TDomain, TAlgebra> {
    public:

        //--------------------------------------------------------------------------------------------------------------

        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;

        using T_IResidualTimeIntegrator = IResidualTimeIntegrator<TDomain, TAlgebra> ;
        using SP_IResidualTimeIntegrator = SmartPtr<T_IResidualTimeIntegrator> ;

        //--------------------------------------------------------------------------------------------------------------

        BasicDriver() : BraidGridFunctionBase<TDomain, TAlgebra>() {}

        BasicDriver(MPI_Comm mpi_temporal, double tstart, double tstop, int steps)
            : BraidGridFunctionBase<TDomain, TAlgebra>(mpi_temporal, tstart, tstop, steps) {
            this->provide_residual = true;
        }

        ~BasicDriver() override = default;

        //--------------------------------------------------------------------------------------------------------------

        int Step(braid_Vector u_, braid_Vector ustop_, braid_Vector fstop_, BraidStepStatus& status) override;;


        int Residual(braid_Vector u_, braid_Vector r_, BraidStepStatus& status) override;

        //--------------------------------------------------------------------------------------------------------------

        void set_approx_space(SmartPtr<ApproximationSpace<TDomain>> spApproxSpace) {
            this->spApproxSpace = spApproxSpace;
        }

        void print_settings() const {}

        void set_default_integrator(SP_IResidualTimeIntegrator integrator) {
            this->default_integrator = integrator;
        }

        void set_integrator(size_t level, SP_IResidualTimeIntegrator integrator);

        SP_IResidualTimeIntegrator get_integrator(size_t level);

        //--------------------------------------------------------------------------------------------------------------

        std::vector<SP_IResidualTimeIntegrator> time_integrator_list;
        SP_IResidualTimeIntegrator default_integrator;

    };



    template<typename TDomain, typename TAlgebra>
    int BasicDriver<TDomain, TAlgebra>::Step(braid_Vector u_, braid_Vector ustop_, braid_Vector fstop_,
        BraidStepStatus &status) {
        int level;
        status.GetLevel(&level);

        int iteration;
        status.GetIter(&iteration);

        double t_start, t_stop;
        status.GetTstartTstop(&t_start, &t_stop);
        std::cout << "times of driver: " << t_start << "\t" << t_stop << std::endl;
        std::cout << std::flush;
        extended_info(std::cout << "BasicDriver::Step --- "
            << " level=" << level
            << " iteration=" << iteration
            << " t_start=" << t_start
            << " t_stop=" << t_stop
            << std::endl << std::flush;)

        auto csp_u_tstop_approx = (*static_cast<SP_GridFunction *>(ustop_->value))->clone();
        auto sp_u_approx_tstart = (*static_cast<SP_GridFunction *>(u_->value))->clone();

        __debug(std::cout << "csp_u_tstop_approx::size = " << csp_u_tstop_approx->size() << std::endl << std::flush);
        __debug(std::cout << "sp_u_approx_tstart::size = " << sp_u_approx_tstart->size() << std::endl << std::flush);

        //[[maybe_unused]] bool success = true;
        auto integrator = this->get_integrator(level);


        if (fstop_ != nullptr) {

            SP_GridFunction fstop = *static_cast<SP_GridFunction *>(fstop_->value);
            integrator->apply(csp_u_tstop_approx, t_stop,
                              sp_u_approx_tstart, t_start,
                              fstop);
        } else {

            bool success = integrator->apply(csp_u_tstop_approx, t_stop,
                                             sp_u_approx_tstart, t_start);
            if (!success) {
                exit(0);
            }

        }

        //extended_info(if(!success) {
        //    this->m_log->o << " [Warning] Convergence was not reached by solver " << std::endl << std::flush;
        //})

        //if (time_adaptivity) {

        //}

        (*static_cast<SP_GridFunction *>(u_->value)) = csp_u_tstop_approx;

        /*{
                int t_index = 0;
                status.GetTIndex(&t_index);
                u_->time = t_stop;
                u_->level_index = t_index+1;
                u_->level = level;
            }*/

        write_script(this->script->Step(u_, ustop_, fstop_, status);)
        return 0;
    }

    template<typename TDomain, typename TAlgebra>
    int BasicDriver<TDomain, TAlgebra>::Residual(braid_Vector u_, braid_Vector r_, BraidStepStatus &status) {
        int level;
        status.GetLevel(&level);

        int iteration;
        status.GetIter(&iteration);

        double t_start, t_stop;
        status.GetTstartTstop(&t_start, &t_stop);


        extended_info( std::cout << "BasicDriver::Residual --- "
            << " level=" << level
            << " iteration=" << iteration
            << " t_start=" << t_start
            << " t_stop=" << t_stop
            << std::endl << std::flush;)

        SP_GridFunction u_tstop = (*static_cast<SP_GridFunction *>(u_->value))->clone();
        SP_GridFunction u_tstart = (*static_cast<SP_GridFunction *>(r_->value))->clone();

        auto integrator = this->get_integrator(level);
        auto result_b = integrator->defect(u_tstop, t_stop, u_tstart, t_start);

        (*result_b) *= -1;
        (*static_cast<SP_GridFunction *>(r_->value)) = result_b;

        /*{
                int t_index = 0;
                status.GetTIndex(&t_index);
                u_->time = t_stop;
                u_->level_index = t_index+1;
                u_->level = level;
            }*/
        write_script(this->script->Residual(u_, r_, status);)

        return 0;
    }

    template<typename TDomain, typename TAlgebra>
    void BasicDriver<TDomain, TAlgebra>::set_integrator(size_t level, SP_IResidualTimeIntegrator integrator) {
        if (this->time_integrator_list.size() < level +1) {
            this->time_integrator_list.resize(level + 1, SPNULL);
        }
        this->time_integrator_list[level] = integrator;
    }

    template<typename TDomain, typename TAlgebra>
    typename BasicDriver<TDomain, TAlgebra>::SP_IResidualTimeIntegrator BasicDriver<TDomain, TAlgebra>::
    get_integrator(size_t level) {
        if (this->time_integrator_list.size() < level +1) {
            return this->default_integrator;
        }
        if (this->time_integrator_list[level] != SPNULL) {
            return this->time_integrator_list[level];
        }
        return this->default_integrator;
    }
}}

#endif