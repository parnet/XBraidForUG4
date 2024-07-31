#ifndef UG_PLUGIN_XBRAIDFORUG4_BASIC_DRIVER_H
#define UG_PLUGIN_XBRAIDFORUG4_BASIC_DRIVER_H

#include "../interface/residual_timeintegrator.h"
#include "gridfunction_base.h"



namespace ug { namespace XBraidForUG4 {

    template <typename TDomain, typename TAlgebra>
    class BasicDriver : public BraidGridFunctionBase<TDomain, TAlgebra> {
    public:

        //--------------------------------------------------------------------------------------------------------------

        typedef GridFunction<TDomain, TAlgebra> T_GridFunction;
        typedef SmartPtr<T_GridFunction> SP_GridFunction;

        typedef IResidualTimeIntegrator<TDomain, TAlgebra> T_IResidualTimeIntegrator;
        typedef SmartPtr<T_IResidualTimeIntegrator> SP_IResidualTimeIntegrator;

        //--------------------------------------------------------------------------------------------------------------

        std::vector<SP_IResidualTimeIntegrator> time_integrator_list;
        SP_IResidualTimeIntegrator default_integrator;

        //--------------------------------------------------------------------------------------------------------------

        BasicDriver() : BraidGridFunctionBase<TDomain, TAlgebra>() {}

        BasicDriver(MPI_Comm mpi_temporal, double tstart, double tstop, int steps)
            : BraidGridFunctionBase<TDomain, TAlgebra>(mpi_temporal, tstart, tstop, steps) {
            this->provide_residual = true;
        }

        ~BasicDriver() override = default;

        //--------------------------------------------------------------------------------------------------------------

        int Step(braid_Vector u_, braid_Vector ustop_, braid_Vector fstop_, BraidStepStatus& pstatus) override {
            int level;
            int idone;
            int tindex;
            int iteration;
            double t_start;
            double t_stop;
            pstatus.GetLevel(&level);
            pstatus.GetTstartTstop(&t_start, &t_stop);
            pstatus.GetDone(&idone);
            pstatus.GetTIndex(&tindex);
            pstatus.GetIter(&iteration);

            std::cout << "BasicDriver::Step --- "
                << " level=" << level
                << " idone=" << idone
                << " tindex=" << tindex
                << " iteration=" << iteration
                << " t_start=" << t_start
                << " t_stop=" << t_stop
                << std::endl << std::flush;

            std::cout << "BasicDriver::Step --- "
                << "Unpack ustop"
                << std::endl << std::flush;
            SP_GridFunction csp_u_tstop_approx = (*(SP_GridFunction*)ustop_->value)->clone();
            std::cout << "BasicDriver::Step --- "
                << "Unpack u_"
                << std::endl << std::flush;
            SP_GridFunction sp_u_approx_tstart = (*(SP_GridFunction*)u_->value)->clone();
            std::cout << "BasicDriver::Step --- "
                << "Get Integrator"
                << std::endl << std::flush;
            bool success = false;
            auto integrator = this->get_integrator(level);


            if (fstop_ != nullptr) {
                std::cout << "BasicDriver::Step --- "
                    << "Unpack fstop"
                    << std::endl << std::flush;
                SP_GridFunction fstop = *(SP_GridFunction*)fstop_->value;
                success = integrator->apply(csp_u_tstop_approx, t_stop,
                                            sp_u_approx_tstart, t_start,
                                            fstop);
            } else {
                std::cout << "BasicDriver::Step --- "
                    << "Apply Integrator"
                    << std::endl << std::flush;
                success = integrator->apply(csp_u_tstop_approx, t_stop,
                                            sp_u_approx_tstart, t_start);
            }
            std::cout << "BasicDriver::Step --- "
                                << "Integrator applied"
                                << std::endl << std::flush;
            if (!success) {
                this->m_log->o << "!!! Failure convergence not reached" << std::endl << std::flush;
            }
            std::cout << "BasicDriver::Step --- "
                    << "Finishes Vector"
                    << std::endl << std::flush;
            (*(SP_GridFunction*)u_->value) = csp_u_tstop_approx;
            std::cout << "BasicDriver::Step --- "
                    << "Vector Set"
                    << std::endl << std::flush;
            return 0;
        };

        int Residual(braid_Vector u_, braid_Vector r_, BraidStepStatus& pstatus) override {
            int level;
            int tindex;
            int iteration;
            double t_start;
            double t_stop;
            pstatus.GetLevel(&level);
            pstatus.GetTIndex(&tindex);
            pstatus.GetIter(&iteration);
            pstatus.GetTstartTstop(&t_start, &t_stop);


            std::cout << "BasicDriver::Residual --- "
                << " level=" << level
                << " tindex=" << tindex
                << " iteration=" << iteration
                << " t_start=" << t_start
                << " t_stop=" << t_stop
                << std::endl << std::flush;

            SP_GridFunction u_tstop = (*(SP_GridFunction*) u_->value)->clone();
            SP_GridFunction u_tstart = (*(SP_GridFunction*) r_->value)->clone();

            auto integrator = this->get_integrator(level);
            auto result_b = integrator->defect(u_tstop, t_stop, u_tstart, t_start);
            (*result_b) *= -1;
            (*(SP_GridFunction*)r_->value) = result_b;
            return 0;
        }

        //--------------------------------------------------------------------------------------------------------------

        void print_settings() {}

        void set_default_integrator(SP_IResidualTimeIntegrator integrator) {
            this->default_integrator = integrator;
        }

        void set_integrator(int level, SP_IResidualTimeIntegrator integrator) {
            if (this->time_integrator_list.size() < level +1) {
                this->time_integrator_list.resize(level + 1, SPNULL);
            }
            this->time_integrator_list[level] = integrator;
        }

        SP_IResidualTimeIntegrator get_integrator(int level) {
            if (this->time_integrator_list.size() < level +1) {
                return this->default_integrator;
            }
            if (this->time_integrator_list[level] != SPNULL) {
                return this->time_integrator_list[level];
            }
            return this->default_integrator;
        }

        //--------------------------------------------------------------------------------------------------------------

    };
}}

#endif //UG_PLUGIN_XBRAIDFORUG4_BASIC_DRIVER_H
