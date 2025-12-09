#ifndef UGPLUGIN_XBRAIDFORUG4_DRIVER_SIMPLE_INTEGRATOR_DRIVER_HPP
#define UGPLUGIN_XBRAIDFORUG4_DRIVER_SIMPLE_INTEGRATOR_DRIVER_HPP


#include <Limex/time_disc/limex_integrator.hpp>

#include "interface/residual_timeintegrator.hpp"
#include "gridfunction_base.hpp"
#include "limex_observer.hpp"


namespace ug{ namespace xbraid {

    template <typename TDomain, typename TAlgebra>
    class SimpleIntegratorDriver final : public BraidGridFunctionBase<TDomain, TAlgebra> {
    public:

        //--------------------------------------------------------------------------------------------------------------

        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;

        using T_LimexObserver = LimexObserver<TDomain, TAlgebra> ;
        using SP_LimexObserver = SmartPtr<T_LimexObserver> ;

        using T_LimexTimeIntegrator = LimexTimeIntegrator<TDomain, TAlgebra> ;
        using SP_LimexTimeIntegrator = SmartPtr<T_LimexTimeIntegrator> ;

        using T_Integrator = SimpleTimeIntegrator<TDomain, TAlgebra>;
        using SP_Integrator = SmartPtr<T_Integrator>;

        using T_Stepper = LinearImplicitEuler<TAlgebra>;

        using T_Solver = IOperatorInverse<typename TAlgebra::vector_type>;
        using SP_Solver = SmartPtr<T_Solver>;

        using T_TimeStepper = LinearImplicitEuler<TAlgebra>;
        using SP_TimeStepper = SmartPtr<T_TimeStepper>;
        //--------------------------------------------------------------------------------------------------------------


        SimpleIntegratorDriver() : BraidGridFunctionBase<TDomain, TAlgebra>() {}


        SimpleIntegratorDriver(MPI_Comm mpi_temporal, double tstart, double tstop, int steps)
            : BraidGridFunctionBase<TDomain, TAlgebra>(mpi_temporal, tstart, tstop, steps) {
            this->provide_residual = false;
        }

        ~SimpleIntegratorDriver() override = default;

        //--------------------------------------------------------------------------------------------------------------

        int Step(braid_Vector u_, braid_Vector ustop_, braid_Vector fstop_, BraidStepStatus& status) override;

        int Residual(braid_Vector u_, braid_Vector r_, BraidStepStatus& status) override;

        int Sync(BraidSyncStatus& status) override;
        //--------------------------------------------------------------------------------------------------------------

        SP_Integrator get_simple_integrator(double dtcurr);


        //--------------------------------------------------------------------------------------------------------------

        /**
         * sets the approximation space for the problem
         * @param sp_approx_space smart pointer of the approximation space
         */
        void set_approx_space(SmartPtr<ApproximationSpace<TDomain>> sp_approx_space) {
            this->sp_approx_space_ = sp_approx_space;
            std::cout << "SimpleIntegratorDriver::used?" << std::endl;
        }


        void print_settings() const {
            std::cout << "SimpleIntegratorDriver::used?" << std::endl;
        }

        void set_integrator(SP_LimexTimeIntegrator integrator);

        void set_tolerance(double loose, double tight);


        number get_level_tolerance(int level) const;

        void set_solver(SP_Solver solver);

        //--------------------------------------------------------------------------------------------------------------

        SP_LimexTimeIntegrator _integrator= SPNULL;
        SP_Integrator _coarse_integrator = SPNULL;
        SP_TimeStepper _time_stepper= SPNULL;
        SP_Solver _solver = SPNULL;

        double _loose = 0.0;
        double _tight = 0.0;

        int _max_level = 2;
        int _current_level = 2;

    };




template<typename TDomain, typename TAlgebra>
int SimpleIntegratorDriver<TDomain, TAlgebra>::Step(braid_Vector u_, braid_Vector ustop_, braid_Vector fstop_,
    BraidStepStatus &status) {

    std::cout << "simple_integrator_driver class called ::step method" << std::endl;

    int level;
    status.GetLevel(&level);

    int iteration;
    status.GetIter(&iteration);

    double t_start, t_stop;
    status.GetTstartTstop(&t_start, &t_stop);
    // double dt = t_stop - t_start;

    double target_tolerance = get_level_tolerance(level);
    std::cout << "set limex target_tolerance = " << target_tolerance << " for level = " << level << std::endl;

    //if (level == 0) // finest level for limex
    //{
    //    _integrator->set_tolerance(target_tolerance);
    //} else { // coarse level for simple integrator
    //}

    int done;
    status.GetDone(&done);


    auto csp_u_tstop_approx = (*static_cast<SP_GridFunction *>(ustop_->value_))->clone();
    auto sp_u_approx_tstart = (*static_cast<SP_GridFunction *>(u_->value_))->clone();

    auto sp_u_approx_tstart_tmp = (*static_cast<SP_GridFunction *>(u_->value_))->clone();

    int index;
    status.GetTIndex(&index);

    //if(done == 1) {

    //} else if (iteration == 0) {
        //SP_LimexObserver observer = make_sp(new T_LimexObserver());
        //_integrator->apply(csp_u_tstop_approx, t_stop, // ø csp_u_tstop_approx -> sp_u_approx_tstart
        //                  sp_u_approx_tstart, t_start);
    //}
    std::cout << "get integrator ----> " << std::endl;
    auto _coarse_integrator = this->get_simple_integrator(t_stop - t_start);
    auto solver = _coarse_integrator->get_solver();
    std::cout << solver->config_string() <<std::endl;
    std::cout << "integrator ready ||--||  " << std::endl;
    _coarse_integrator->apply(csp_u_tstop_approx, t_stop, sp_u_approx_tstart, t_start);
    std::cout << "integrator after apply #---#  " << std::endl;
    //size_t steps = _integrator->get_step() - 1;
    // notify_finalize_step( u, limex step, t, dt)

    (*static_cast<SP_GridFunction *>(u_->value_)) = csp_u_tstop_approx;

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    std::cout << "x_step_output" << std::endl;
    std::cout << "iteration = " << iteration <<std::endl;
    std::cout << "level = " << level <<std::endl;
    std::cout << "tolerance = " << target_tolerance <<std::endl;
    std::cout << "t_index = " << index <<std::endl;
    //std::cout << "steps = " << steps <<std::endl;

    if(level == 0) {
        // int r_factor = static_cast<int>(steps) / 2;
        // if (r_factor < 0) {
        //    r_factor = 1;
        //}
        //status.SetRFactor(r_factor);
    }

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    write_script(this->script_->Step(u_, ustop_, fstop_, status);)
    return 0;
}


template<typename TDomain, typename TAlgebra>
int SimpleIntegratorDriver<TDomain, TAlgebra>::Residual(braid_Vector u_, braid_Vector r_, BraidStepStatus &status) {
    std::cout << "SimpleIntegratorDriver::Residual is not supported. check configuration" << std::endl;
    exit(1);
}

template<typename TDomain, typename TAlgebra>
int SimpleIntegratorDriver<TDomain, TAlgebra>::Sync(BraidSyncStatus& status) {
        std::cout << "SimpleIntegratorDriver::used?" << std::endl;
        __debug(std::cout << "SimpleIntegratorDriver::Sync" << std::endl);
        this->iteration_ += 1;
        write_script(this->script_->Sync(status);)

        return 0;
    }

template<typename TDomain, typename TAlgebra>
typename SimpleIntegratorDriver<TDomain, TAlgebra>::SP_Integrator SimpleIntegratorDriver<TDomain, TAlgebra>::
get_simple_integrator(double dtcurr) {

    if (_coarse_integrator == SPNULL) {
        SP_TimeStepper stepper = _integrator->get_time_stepper(0);
        if (stepper == SPNULL) {
            std::cout << "stepper is nullptr" << std::endl;
        }
        _coarse_integrator = SmartPtr<T_Integrator>(new T_Integrator(stepper));

        //SP_Solver solver = _integrator->get_solver(0);
        //if (solver == SPNULL) {
        //    std::cout << "solver is nullptr" << std::endl;
        //}
        _coarse_integrator->set_solver(this->_solver);

        if (_coarse_integrator == SPNULL) {
            std::cout << "integrator is nullptr" << std::endl;
        }
        /*integrator.set_dt_min(dtcurr/m_vSteps[i]);
                integrator.set_dt_max(dtcurr/m_vSteps[i]);*/
        _coarse_integrator->set_reduction_factor(0.0);                 // quit immediately, if step fails

    }

    SP_GridFunction derivative = this->_integrator->get_time_derivative();
    if (derivative == SPNULL) {
        std::cout << "derivative is nullptr" << std::endl;
    }
    _coarse_integrator->set_derivative(derivative);
    _coarse_integrator->set_time_step(dtcurr);
    _coarse_integrator->set_dt_min(dtcurr); // /(log(m_epsmin)/log(m_tol))
    _coarse_integrator->set_dt_max(dtcurr); // *log(m_epsmin)/log(m_tol)

    auto banach_space = _integrator->get_space();
    if (banach_space == SPNULL) {
        std::cout << "banach_space is nullptr" << std::endl;
    }
    _coarse_integrator->set_banach_space(banach_space);
    return _coarse_integrator;
};



template<typename TDomain, typename TAlgebra>
void SimpleIntegratorDriver<TDomain, TAlgebra>::set_integrator(SP_LimexTimeIntegrator integrator) {
     this->_integrator = integrator;
    }

template<typename TDomain, typename TAlgebra>
void SimpleIntegratorDriver<TDomain, TAlgebra>::set_tolerance(double loose, double tight) {
    std::cout << "loose=" << loose << "is ignored"<< std::endl;
    std::cout << "tight=" << tight << std::endl;
    std::cout << std::endl;
    this->_loose = loose;
    this->_tight = tight;
}

template<typename TDomain, typename TAlgebra>
number SimpleIntegratorDriver<TDomain, TAlgebra>::get_level_tolerance(int level) const {
    int fine_level = 0;
    int base_level = 2;
    double log_loose = log(_loose);
    double log_tight = log(_tight);
    int number_of_level = base_level - fine_level + 1;
    double linear_ratio = static_cast<double>(base_level - level) / static_cast<double>( number_of_level -1 );
    double linear_interpolate = log_loose + linear_ratio * (log_tight - log_loose);
    return exp(linear_interpolate);
};

template<typename TDomain, typename TAlgebra>
void SimpleIntegratorDriver<TDomain, TAlgebra>::set_solver(SP_Solver solver) {
        this->_solver = solver;
    };



}}

#endif
