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

        using T_DebugWriter = IDebugWriter<TAlgebra>;
        using SP_DebugWriter = SmartPtr<T_DebugWriter>;
        //--------------------------------------------------------------------------------------------------------------


        SimpleIntegratorDriver() : BraidGridFunctionBase<TDomain, TAlgebra>() {
        }


        SimpleIntegratorDriver(MPI_Comm mpi_temporal, double tstart, double tstop, int steps)
            : BraidGridFunctionBase<TDomain, TAlgebra>(mpi_temporal, tstart, tstop, steps) {
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

        void set_debug_write(SP_DebugWriter debug_writer ) {
            this->_debug_writer = debug_writer;
        }
        //--------------------------------------------------------------------------------------------------------------

        SP_LimexTimeIntegrator _limex_integrator= SPNULL;
        SP_Integrator _simple_integrator = SPNULL;
        SP_TimeStepper _time_stepper= SPNULL;
        SP_Solver _solver = SPNULL;
        SP_DebugWriter _debug_writer = SPNULL;

        double _loose = 0.0;
        double _tight = 0.0;

        int _max_level = 2;
        int _current_level = 2;
        int _gridstep = 2;
    };




template<typename TDomain, typename TAlgebra>
int SimpleIntegratorDriver<TDomain, TAlgebra>::Step(braid_Vector u_, braid_Vector ustop_,
    braid_Vector fstop_, BraidStepStatus &status) {
    int level;
    status.GetLevel(&level);
    //std::cout << "SimpleIntegratorDriver::Step " << level << std::endl;
    double t_start, t_stop;
    status.GetTstartTstop(&t_start, &t_stop);
    //auto csp_u_tstop_approx =(*(SP_GridFunction*)(ustop_->value_))->clone();
    auto csp_u_tstop_approx =(*static_cast<SP_GridFunction *>(u_->value_))->clone();
    auto sp_u_approx_tstart = (*static_cast<SP_GridFunction *>(u_->value_))->clone();

    //std::cout << "u_n : " << csp_u_tstop_approx.get() << " norm="<< this->norm_->norm(csp_u_tstop_approx) << std::endl;
    //std::cout << "u_0 : " << sp_u_approx_tstart.get() << " norm="<< this->norm_->norm(sp_u_approx_tstart) << std::endl;

    double dt = (t_stop - t_start);
    //std::cout << "Integrating from: " << t_start << " to: " << t_stop << " with dt = " << dt << std::endl;
    auto simple_integrator = this->get_simple_integrator(dt);
    //simple_integrator->set_debug(_debug_writer);
    simple_integrator->apply(csp_u_tstop_approx, t_stop,
                             sp_u_approx_tstart, t_start);
    (*static_cast<SP_GridFunction *>(u_->value_)) = csp_u_tstop_approx;
    //((SP_GridFunction*)(u_->value_))->operator=( csp_u_tstop_approx);
    //std::cout << "r : " << (*(SP_GridFunction*)(u_->value_)).get() << " norm: " << this->norm_->norm(csp_u_tstop_approx)  << std::endl;
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
    if (_simple_integrator == SPNULL) {
        SP_TimeStepper stepper = _limex_integrator->get_time_stepper(0);
        if (stepper == SPNULL) {
            std::cout << "stepper is nullptr" << std::endl;
        }
        auto banach_space = _limex_integrator->get_space();
        if (banach_space == SPNULL) {
            std::cout << "banach_space is nullptr" << std::endl;
        }
        stepper->set_matrix_cache(false);
        _simple_integrator = SmartPtr<T_Integrator>(new T_Integrator(stepper,banach_space));
        auto solver = _limex_integrator->get_solver(0);
        _simple_integrator->set_solver(solver);
        if (_simple_integrator == SPNULL) {
            std::cout << "integrator is nullptr" << std::endl;
        }
        _simple_integrator->set_reduction_factor(0.0);
    }
    _simple_integrator->set_time_step(dtcurr);
    _simple_integrator->set_dt_min(dtcurr);
    _simple_integrator->set_dt_max(dtcurr);
    return _simple_integrator;
};



template<typename TDomain, typename TAlgebra>
void SimpleIntegratorDriver<TDomain, TAlgebra>::set_integrator(SP_LimexTimeIntegrator integrator) {
     this->_limex_integrator = integrator;
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
        std::cout << "recv - config string : " << std::endl;
        std::cout << this->_solver->config_string() << std::endl;
    };



}}

#endif
