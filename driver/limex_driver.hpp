#ifndef UGPLUGIN_XBRAIDFORUG4_DRIVER_LIMEX_DRIVER_HPP
#define UGPLUGIN_XBRAIDFORUG4_DRIVER_LIMEX_DRIVER_HPP


#include <boost/graph/graph_traits.hpp>
#include <Limex/time_disc/limex_integrator.hpp>

#include "interface/residual_timeintegrator.hpp"
#include "gridfunction_base.hpp"
#include "limex_observer.hpp"


namespace ug{ namespace xbraid {

    template <typename TDomain, typename TAlgebra>
    class LimexDriver final : public BraidGridFunctionBase<TDomain, TAlgebra> {
    public:

        //--------------------------------------------------------------------------------------------------------------

        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;

        using T_LimexObserver = LimexObserver<TDomain, TAlgebra> ;
        using SP_LimexObserver = SmartPtr<T_LimexObserver> ;

        using T_LimexTimeIntegrator = LimexTimeIntegrator<TDomain, TAlgebra> ;
        using SP_LimexTimeIntegrator = SmartPtr<T_LimexTimeIntegrator> ;

        //--------------------------------------------------------------------------------------------------------------


        LimexDriver() : BraidGridFunctionBase<TDomain, TAlgebra>() {}


        LimexDriver(MPI_Comm mpi_temporal, double tstart, double tstop, int steps)
            : BraidGridFunctionBase<TDomain, TAlgebra>(mpi_temporal, tstart, tstop, steps) {
            this->provide_residual = false;
        }

        ~LimexDriver() override = default;

        //--------------------------------------------------------------------------------------------------------------

        int Step(braid_Vector u_, braid_Vector ustop_, braid_Vector fstop_, BraidStepStatus& status) override;

        int Residual(braid_Vector u_, braid_Vector r_, BraidStepStatus& status) override;

        int Sync(BraidSyncStatus& status) override;
        //--------------------------------------------------------------------------------------------------------------

        /**
         * sets the approximation space for the problem
         * @param sp_approx_space smart pointer of the approximation space
         */
        void set_approx_space(SmartPtr<ApproximationSpace<TDomain>> sp_approx_space) {
            this->sp_approx_space_ = sp_approx_space;
            std::cout << "LimexDriver::used?" << std::endl;
        }


        void print_settings() const {
            std::cout << "LimexDriver::used?" << std::endl;
        }

        void set_integrator(SP_LimexTimeIntegrator integrator);

        void set_tolerance(double loose, double tight) {
            std::cout << "loose=" << loose << std::endl;
            std::cout << "tight=" << tight << std::endl;
            std::cout << std::endl;
            this->_loose = loose;
            this->_tight = tight;
        }


        number get_level_tolerance(int level) {
            int fine_level = 0;
            int base_level = 2; // todo move
            double log_loose = log(_loose);
            double log_tight = log(_tight);
            int number_of_level = base_level - fine_level + 1;
            double linear_ratio = static_cast<double>(base_level - level) / static_cast<double>( number_of_level -1 );
            double linear_interpolate = log_loose + linear_ratio * (log_tight - log_loose);
            return exp(linear_interpolate);
        }
        //--------------------------------------------------------------------------------------------------------------

        SP_LimexTimeIntegrator _integrator;
        double _loose = 0.0;
        double _tight = 0.0;

        int _max_level = 2;
        int _current_level = 2;

    };




template<typename TDomain, typename TAlgebra>
int LimexDriver<TDomain, TAlgebra>::Step(braid_Vector u_, braid_Vector ustop_, braid_Vector fstop_,
    BraidStepStatus &status) {

    int level;
    status.GetLevel(&level);

    int iteration;
    status.GetIter(&iteration);

    double t_start, t_stop;
    status.GetTstartTstop(&t_start, &t_stop);
    double dt = t_stop - t_start;

    double target_tolerance = get_level_tolerance(level);
    std::cout << "set limex target_tolerance = " << target_tolerance << " for level = " << level << std::endl;

    _integrator->set_tolerance(target_tolerance);
    if (level > 0 ){
        _integrator->set_time_step(dt/2);
    }
    int done;
    status.GetDone(&done);

    auto core = status.GetCore();


    auto csp_u_tstop_approx = (*static_cast<SP_GridFunction *>(ustop_->value_))->clone();
    auto sp_u_approx_tstart = (*static_cast<SP_GridFunction *>(u_->value_))->clone();

    auto sp_u_approx_tstart_tmp = (*static_cast<SP_GridFunction *>(u_->value_))->clone();

    int index;
    int index_mod;
    status.GetTIndex(&index);

    //index_mod = 100'000 * (iteration+1) + 10'000*(level+1) + index;
    //this->out_->step_process(sp_u_approx_tstart_tmp,index_mod,t_stop, t_stop-t_start);

    //double value = this->norm_->norm(sp_u_approx_tstart_tmp);
    //std::cout << "norm before apply: " << value << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    // todo prepare integration
    if(done == 1) {
        // todo attach output observer
    } else if (iteration == 0) {
        //SP_LimexObserver observer = make_sp(new T_LimexObserver());
    }


    _integrator->apply(csp_u_tstop_approx, t_stop, // ø csp_u_tstop_approx -> sp_u_approx_tstart
                      sp_u_approx_tstart, t_start);

    size_t steps = _integrator->get_step() - 1;
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
    std::cout << "steps = " << steps <<std::endl;

    if(level == 0) {
        int r_factor = steps / 2;
        if (r_factor < 0) {
            r_factor = 1;
        }
        status.SetRFactor(r_factor);
    }

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    write_script(this->script_->Step(u_, ustop_, fstop_, status);)
    return 0;
}


template<typename TDomain, typename TAlgebra>
int LimexDriver<TDomain, TAlgebra>::Residual(braid_Vector u_, braid_Vector r_, BraidStepStatus &status) {
    std::cout << "LimexDriver::Residual is not supported. check configuration" << std::endl;
    exit(1);
    return 0;
}

template<typename TDomain, typename TAlgebra>
int LimexDriver<TDomain, TAlgebra>::Sync(BraidSyncStatus& status) {
        std::cout << "LimexDriver::used?" << std::endl;
        __debug(std::cout << "LimexDriver::Sync" << std::endl);
        this->iteration_ += 1;
        write_script(this->script_->Sync(status);)

        return 0;
    };



template<typename TDomain, typename TAlgebra>
void LimexDriver<TDomain, TAlgebra>::set_integrator(SP_LimexTimeIntegrator integrator) {
     this->_integrator = integrator;
    };




}}

#endif
