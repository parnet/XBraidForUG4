#ifndef UGPLUGIN_XBRAIDFORUG4_DUMMY_SIMPLE_INTEGRATOR_DRIVER_HPP
#define UGPLUGIN_XBRAIDFORUG4_DUMMY_SIMPLE_INTEGRATOR_DRIVER_HPP


#include <Limex/time_disc/limex_integrator.hpp>

#include "interface/residual_timeintegrator.hpp"
#include "../driver/gridfunction_base.hpp"
#include "../driver/limex_observer.hpp"

/*
namespace ug{ namespace xbraid {

template <typename TDomain, typename TAlgebra>
class IsolatedSimpleIntegratorDriver final  {
public:
        /*using SP_SpaceTimeCommunicator = SmartPtr<SpaceTimeCommunicator> ;

        using T_DomainDisc = IDomainDiscretization<TAlgebra> ;
        using SP_DomainDisc = SmartPtr<T_DomainDisc> ;

        using T_VectorValueType = typename TAlgebra::vector_type::value_type;
        * /
        using T_GridFunction= GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;

        /*using T_ITimeIntegrator = ITimeIntegrator<TDomain, TAlgebra> ;
        using SP_TimeIntegrator = SmartPtr<T_ITimeIntegrator> ;*/

        /*using T_BraidInitializer = BraidInitializer<TDomain, TAlgebra> ;
        using SP_BraidInitializer = SmartPtr<T_BraidInitializer> ;*/

        /*using T_ITimeIntegratorObserver = ITimeIntegratorObserver<TDomain, TAlgebra> ;
        using SP_IObserver = SmartPtr<T_ITimeIntegratorObserver> ;

        using T_IXBraidTimeIntegratorObserver = IXBraidTimeIntegratorObserver<TDomain, TAlgebra> ;
        using SP_IXBraidTimeIntegratorObserver = SmartPtr<T_IXBraidTimeIntegratorObserver> ;

        using T_ParallelLogger = ParallelLogger ;
        using SP_ParallelLogger = SmartPtr<T_ParallelLogger> ;

        using T_ApproximationSpace = ApproximationSpace<TDomain>;
        using SP_ApproximationSpace = SmartPtr<T_ApproximationSpace> ;

        using T_SpatialNorm = BraidSpatialNorm<TDomain, TAlgebra> ;
        using SP_SpatialNorm = SmartPtr<T_SpatialNorm> ;

        using T_BraidWriteScript = BraidWriteScript<TDomain, TAlgebra> ;
        using SP_BraidWriteScript = SmartPtr<T_BraidWriteScript> ;


        using T_LimexObserver = LimexObserver<TDomain, TAlgebra> ;
        using SP_LimexObserver = SmartPtr<T_LimexObserver> ;* /

        using T_LimexTimeIntegrator = LimexTimeIntegrator<TDomain, TAlgebra> ;
        using SP_LimexTimeIntegrator = SmartPtr<T_LimexTimeIntegrator> ;

        using T_Integrator = SimpleTimeIntegrator<TDomain, TAlgebra>;
        using SP_Integrator = SmartPtr<T_Integrator>;

        /*using T_Stepper = LinearImplicitEuler<TAlgebra>;
        using SP_Stepper = SmartPtr<T_Stepper>;

        using T_Solver = IOperatorInverse<typename TAlgebra::vector_type>;
        using SP_Solver = SmartPtr<T_Solver>;* /

        using T_TimeStepper = LinearImplicitEuler<TAlgebra>;
        using SP_TimeStepper = SmartPtr<T_TimeStepper>;

        using T_DebugWriter = IDebugWriter<TAlgebra>;
        using SP_DebugWriter = SmartPtr<T_DebugWriter>;

        /*using T_SpatialGridTransfer = SpatialGridTransfer<TDomain,TAlgebra>;
        using SP_SpatialGridTransfer = SmartPtr<T_SpatialGridTransfer>;* /


    IsolatedSimpleIntegratorDriver() = default;

    /*IsolatedSimpleIntegratorDriver(double tstart, double tstop, int steps) {
        this->_tstart = tstart;
        this->_tstop = tstop;
        this->_steps = steps;
        }* /


    ~IsolatedSimpleIntegratorDriver() = default;

    auto step(SP_GridFunction u_, SP_GridFunction ustop_, double t_start, double t_stop) -> SP_GridFunction;

    /*int Residual(braid_Vector u_, braid_Vector r_, BraidStepStatus& status);* /


    SP_Integrator get_simple_integrator(double dtcurr);

    /*void set_approx_space(SP_ApproximationSpace sp_approx_space) {
            this->sp_approx_space_ = sp_approx_space;
        }* /

    void set_integrator(SP_LimexTimeIntegrator integrator);

    /*void set_tolerance(double loose, double tight);*/


    /*number get_level_tolerance(int level) const;*/

    /*void set_solver(SP_Solver solver);* /

    void set_debug_write(SP_DebugWriter debug_writer ) {
        this->_debug_writer = debug_writer;
    }

    /*int Init(double t, SP_GridFunction u_ptr);*/

    /*int Clone(braid_Vector u_, braid_Vector* v_ptr);*/

    /*int Free(braid_Vector u_);*/

    // y = alpha * x + beta*y
    /*int Sum(double alpha, braid_Vector x_, double beta, braid_Vector y_);*/

    /*int SpatialNorm(braid_Vector u_, double* norm_ptr);*/

    /*int Access(braid_Vector u_, BraidAccessStatus& status);*/

    /*int BufSize(int* size_ptr, BraidBufferStatus& status);

    int BufPack(braid_Vector u_, void* buffer, BraidBufferStatus& status);

    int BufUnpack(void* buffer, braid_Vector* u_ptr, BraidBufferStatus& status);

    void pack(void* buffer, T_GridFunction* u_ref, int* buffer_size);

    void unpack(void* buffer, T_GridFunction* u_ref, int* buffer_size);*/

    /*int Sync(BraidSyncStatus& status);*/

    /*void set_spatial_grid_transfer(SP_SpatialGridTransfer spatial_grid_transfer) {
        this->spatial_grid_transfer = spatial_grid_transfer;
    }

    int Coarsen(braid_Vector fu_,
              braid_Vector* cu_ptr,
              BraidCoarsenRefStatus &status);

    int Refine(braid_Vector cu_, braid_Vector *fu_ptr, BraidCoarsenRefStatus &status);*/

    /*void set_paralog(SP_ParallelLogger log) {
        this->log_ = log;
    }*/

    /*void init() {
        __send_recv_times(this->timer = BraidTimer(););
        this->log_->init();
        write_script(this->script_ = make_sp(new T_BraidWriteScript(this->comm_));) // å
    }*/

    /*void set_time_values(double startTime, double endTime, int n) {
        this->_tstart = startTime;
        this->_tstop = endTime;
        this->_steps = n;
    }*/

    /*void set_start_time(double startTime) {
        this->_tstart = startTime;
    }*/

    /*void set_end_time(double endTime) {
        this->_tstop = endTime;
    }*/

    /*void set_number_of_timesteps(int n) {
        this->_steps = n;
    }*/

    /*void set_start_vector(SP_GridFunction p_u0) {
        this->u0_ = p_u0;
    }*/

    /*void attach_xbraid_observer(SP_IXBraidTimeIntegratorObserver p_out) {
        this->xb_out_ = p_out;
    }*/

    /*void attach_observer(SP_IObserver p_out) {
        this->out_ = p_out;
    }*/

    /*void set_norm_provider(SP_SpatialNorm norm) {
        this->norm_ = norm;
    }*/

    /*void set_max_levels(size_t levelcount) {
        this->levels_ = levelcount;
    }*/

    /*void set_initializer(SP_BraidInitializer initializer) {
        std::cout << "set-initializer" << std::endl;
        this->initializer_ = initializer;
    }*/

    /*void set_domain(SP_DomainDisc domain) {
        this->domain_disc_ = domain;
    }*/

/*#ifdef FEATURE_SPATIAL_REFINE
    void set_level_num_ref(size_t level, int num_ref);
#endif* /

        //size_t m_init_counter = 0;

        //SP_DomainDisc domain_disc_;
        //SP_SpaceTimeCommunicator comm_;
        //SP_ParallelLogger log_;
        //SP_GridFunction u0_; // used for buffer and for rhs (residual) construction,

        //SP_IObserver out_;
        //SP_IXBraidTimeIntegratorObserver xb_out_;
        //SP_BraidInitializer initializer_;
        //SP_SpatialNorm norm_;
        //bool can_residual_method_ = false;// error estimation and refine
        //bool refine_time_ = false;// error estimation and refine
        //bool restimate_ = false;
        //int levels_ = 0;
        //int iteration_ = 0;

        //int norm_counter = 0;
        //PIOGridFunction<TDomain,TAlgebra> pio_grid_function_;

        SP_LimexTimeIntegrator _limex_integrator = SPNULL;
        SP_Integrator _simple_integrator = SPNULL;
        //SP_TimeStepper _time_stepper= SPNULL;
        //SP_Solver _solver = SPNULL;
        SP_DebugWriter _debug_writer = SPNULL;

        //double _tstart;
        //double _tstop;
        //double _steps;

        //double _loose = 0.0;
        //double _tight = 0.0;

        //int _max_level = 2;
        //int _current_level = 2;
        //int _gridstep = 2;

        //__send_recv_times(BraidTimer timer_;)
        //write_script(SP_BraidWriteScript script_;)

#ifdef FEATURE_SPATIAL_REFINE
        //std::vector<int> level_num_ref;
        //SP_ApproximationSpace sp_approx_space_ = SPNULL;
        //SP_SpatialGridTransfer spatial_grid_transfer;
#endif
    };


// ---------------------------------------------------------------------------------------------------------------------

template<typename TDomain, typename TAlgebra>
auto IsolatedSimpleIntegratorDriver<TDomain, TAlgebra>::
step(SP_GridFunction u_, SP_GridFunction ustop_, double t_start, double t_stop ) -> SP_GridFunction {
    auto csp_u_tstop_approx = ustop_->clone();
    auto sp_u_approx_tstart = u_->clone();
    double dt = (t_stop - t_start);
    std::cout << "Integrating from: " << t_start << " to: " << t_stop << " with dt = " << dt << std::endl;
    auto _coarse_integrator = this->get_simple_integrator(t_stop - t_start);
    _coarse_integrator->set_debug(_debug_writer);
    _coarse_integrator->apply(csp_u_tstop_approx, t_stop,
                              sp_u_approx_tstart, t_start);
    //}


    return csp_u_tstop_approx;
}

/*
template<typename TDomain, typename TAlgebra>
int IsolatedSimpleIntegratorDriver<TDomain, TAlgebra>::
Residual(braid_Vector u_, braid_Vector r_, BraidStepStatus &status) {
    std::cout << "IsolatedSimpleIntegratorDriver::Residual is not supported. check configuration" << std::endl;
    exit(1);
}*/

/*template<typename TDomain, typename TAlgebra>
int IsolatedSimpleIntegratorDriver<TDomain, TAlgebra>::
Sync(BraidSyncStatus& status) {
        std::cout << "IsolatedSimpleIntegratorDriver::used?" << std::endl;
        __debug(std::cout << "IsolatedSimpleIntegratorDriver::Sync" << std::endl);
        this->iteration_ += 1;
        write_script(this->script_->Sync(status);)

        return 0;
}* /


template<typename TDomain, typename TAlgebra>
typename IsolatedSimpleIntegratorDriver<TDomain, TAlgebra>::SP_Integrator IsolatedSimpleIntegratorDriver<TDomain, TAlgebra>::
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
void IsolatedSimpleIntegratorDriver<TDomain, TAlgebra>::set_integrator(SP_LimexTimeIntegrator integrator) {
     this->_limex_integrator = integrator;
    }

/*template<typename TDomain, typename TAlgebra>
void IsolatedSimpleIntegratorDriver<TDomain, TAlgebra>::set_tolerance(double loose, double tight) {
    std::cout << "loose=" << loose << "is ignored"<< std::endl;
    std::cout << "tight=" << tight << std::endl;
    std::cout << std::endl;
    this->_loose = loose;
    this->_tight = tight;
}*/

/*template<typename TDomain, typename TAlgebra>
number IsolatedSimpleIntegratorDriver<TDomain, TAlgebra>::get_level_tolerance(int level) const {
    int fine_level = 0;
    int base_level = 2;
    double log_loose = log(_loose);
    double log_tight = log(_tight);
    int number_of_level = base_level - fine_level + 1;
    double linear_ratio = static_cast<double>(base_level - level) / static_cast<double>( number_of_level -1 );
    double linear_interpolate = log_loose + linear_ratio * (log_tight - log_loose);
    return exp(linear_interpolate);
};*/

/*template<typename TDomain, typename TAlgebra>
void IsolatedSimpleIntegratorDriver<TDomain, TAlgebra>::set_solver(SP_Solver solver) {
        this->_solver = solver;
        std::cout << "recv - config string : " << std::endl;
        std::cout << this->_solver->config_string() << std::endl;
    };*/


/*
template<typename TDomain, typename TAlgebra>
int IsolatedSimpleIntegratorDriver<TDomain, TAlgebra>::Refine(braid_Vector cu_,
                                     braid_Vector *fu_ptr,
                                     BraidCoarsenRefStatus &status)
        {
#ifdef FEATURE_SPATIAL_REFINE
            int mgrit_level_fine;
            status.GetLevel(&mgrit_level_fine);
            int mgrit_level_coarse = mgrit_level_fine+1;
            std::cout << "MGRIT - Refining: " <<  mgrit_level_coarse << " ---> " << mgrit_level_fine << "    ---    "<< std::endl <<std::flush;
            __debug(std::cout << "MGRIT - Refining: " <<  mgrit_level_coarse << " ---> " << mgrit_level_fine << "    ---    "<< std::endl <<std::flush);

            const int gmg_level_fine = level_num_ref[mgrit_level_fine];
            const int gmg_level_coarse = level_num_ref[mgrit_level_coarse];
            __debug(std::cout << "Grid NumRef: " <<  gmg_level_coarse<< " ---> " <<  gmg_level_fine<< "    ---    " << std::endl <<std::flush);

            if ( gmg_level_fine == gmg_level_coarse){ // no refinement
                __debug(std::cout <<  "not refine --> clone " << std::endl);
                this->Clone(cu_, fu_ptr);
                __debug(std::cout <<  "not refining --> Finished " << std::endl <<std::flush);
            } else {
                SP_GridFunction sp_cu = (*static_cast<SP_GridFunction *>((cu_)->value_));

                const size_t refs = gmg_level_fine - gmg_level_coarse;


                //std::stringstream filename;
                //filename << "it_"<< this->iteration << "_refine_before_u_" << cu_->index;
                //pio_grid_function_.write(sp_cu,filename.str().c_str());


                SP_GridFunction tmp = sp_cu;
                for ( size_t i = 0; i < refs; ++i) {
                    std::cout << "refining step: " << i << std::flush << std::endl;
                    tmp = this->spatial_grid_transfer->prolongate(tmp);
                }
                std::cout <<  " ----------------------  =" << gmg_level_fine << std::endl <<std::flush;

                auto * sp_fu = new SmartPtr<T_GridFunction>(tmp);
                //sp_fu->get()->set_storage_type(sp_cu->get_storage_mask());
                //sp_fu->enable_redistribution(sp_cu->redistribution_enabled());

                __debug(std::cout <<  "prolongation: done! "<< std::endl <<std::flush);

                auto* u = static_cast<BraidVector *>(malloc(sizeof(BraidVector)));
                u->value_ = sp_fu;
                *fu_ptr = u;
                // (*fu_ptr)->index = indexpool++;

                //std::stringstream filename_after;
                //filename_after << "it_"<< this->iteration << "_refine_after_u_" << (*fu_ptr)->index;
                //pio_grid_function_.write(*sp_fu,filename_after.str().c_str());
            }
            / *{
                int t_index;
                status.GetTIndex(&t_index);
                (*fu_ptr)->time = cu_->time;
                (*fu_ptr)->level_index = t_index;
                (*fu_ptr)->level = cu_->level+1;
            }* /
            write_script(this->script_->Refine(cu_, fu_ptr, status);)
            return 0;
#else
            this->Clone(cu_, fu_ptr); // no refinement
            //void (* function)(T_GridFunction &, const T_GridFunction &));
            write_script(this->script->Coarsen(cu_, fu_ptr, status);)
            return 0;
#endif

        }


template<typename TDomain, typename TAlgebra>
int IsolatedSimpleIntegratorDriver<TDomain, TAlgebra>::Coarsen(braid_Vector fu_,
                          braid_Vector* cu_ptr,
                          BraidCoarsenRefStatus &status)
        {
#ifdef FEATURE_SPATIAL_REFINE
            int mgrit_level_fine;
            status.GetLevel(&mgrit_level_fine);

            const int mgrit_level_coarse = mgrit_level_fine+1;
            std::cout << "MGRIT - Coarsening: " << mgrit_level_fine  << " ---> " << mgrit_level_coarse  << "    ---    "<< std::endl <<std::flush
            __debug(std::cout << "MGRIT - Coarsening: " <<mgrit_level_fine  << " ---> " <<  mgrit_level_coarse  << "    ---    "<< std::endl <<std::flush);
            // ---------------------------------------------------------------------------------------------------------


            int gmg_level_fine = level_num_ref[mgrit_level_fine];
            int gmg_level_coarse = level_num_ref[mgrit_level_coarse];
            __debug(std::cout << "Grid NumRef: " << gmg_level_fine << " ---> " << gmg_level_coarse<< "    ---    " << std::endl <<std::flush);
            if ( gmg_level_fine == gmg_level_coarse){ // no coarsening
                __debug(std::cout <<  "not coarsen --> clone " << std::endl <<std::flush);
                this->Clone(fu_, cu_ptr);
                __debug(std::cout <<  "not coarsen --> Finished " << std::endl <<std::flush);
            } else {
                SP_GridFunction sp_fu = (*static_cast<SP_GridFunction *>((fu_)->value_));
                const size_t refs = gmg_level_fine - gmg_level_coarse;

                //std::stringstream filename;
                //filename << "it_"<< this->iteration << "_coarsen_before_u_" << fu_->index;
                //pio_grid_function_.write(sp_fu,filename.str().c_str());


                SP_GridFunction tmp = sp_fu;

                for ( size_t i = 0; i < refs; ++i) {
                    std::cout << "coarsen step: " << i << std::flush << std::endl;
                    tmp = this->spatial_grid_transfer->restrict(tmp);
                    auto result = tmp->clone();

                }

                std::cout <<  " ----------------------  =" << gmg_level_coarse << std::endl <<std::flush;
                auto * sp_cu = new SmartPtr<T_GridFunction>(tmp);
                //sp_cu->get()->set_storage_type(sp_fu->get_storage_mask());
                //sp_cu->enable_redistribution(sp_fu->redistribution_enabled());
                __debug(std::cout <<  "restriction: done! "<< std::endl <<std::flush);

                auto* u = (BraidVector*)malloc(sizeof(BraidVector));
                u->value_ = sp_cu;
                *cu_ptr = u;
                //(*cu_ptr)->index = indexpool++;

                //std::stringstream filename_after;
                //filename_after << "it_"<< this->iteration << "_coarsen_after_u_" << (*cu_ptr)->index;
                //pio_grid_function_.write(*sp_cu,filename_after.str().c_str());

            }
            / *{
                int t_index;
                status.GetTIndex(&t_index);
                (*cu_ptr)->time = fu_->time;
                (*cu_ptr)->level_index = t_index;
                (*cu_ptr)->level = fu_->level+1;
            }* /
            write_script(this->script_->Coarsen(fu_, cu_ptr, status);)
            __debug(std::cout << "~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ " << std::endl << std::flush);
            return 0;
#else
            this->Clone(fu_, cu_ptr); // no coarsening
            write_script(this->script->Coarsen(fu_, cu_ptr, status);)
            return 0;
#endif
        }
*/
/*
template<typename TDomain, typename TAlgebra>
int IsolatedSimpleIntegratorDriver<TDomain, TAlgebra>::Clone(braid_Vector u_, braid_Vector* v_ptr) {
    __debug(std::cout << "GridFunctionBaseDriver::Clone" << std::endl);

    auto* v = static_cast<BraidVector *>(malloc(sizeof(BraidVector)));
    auto* uref = static_cast<SP_GridFunction *>(u_->value_);
    auto* vref = new SP_GridFunction();
    *vref = uref->get()->clone();
    v->value_ = vref;
    v->time_ = u_->time_;
    *v_ptr = v;

    / *{
        v->time = u_->time;
        v->level_index = u_->level_index;
        v->level = u_->level;
    }* /
    write_script(this->script_->Clone(u_,v_ptr);)
    return 0;
};*/


/*template<typename TDomain, typename TAlgebra>
int IsolatedSimpleIntegratorDriver<TDomain, TAlgebra>::SpatialNorm(braid_Vector u_, double* norm_ptr) {
        __debug(std::cout << "GridFunctionBaseDriver::SpatialNorm" << std::endl);
        *norm_ptr = 0;
        auto* uref = static_cast<SP_GridFunction *>(u_->value_);
        SP_GridFunction tempobject_output = uref->get()->clone();

        SP_GridFunction tempobject = uref->get()->clone();
        *norm_ptr = norm_->norm(tempobject);

        //#ifdef FEATURE_WRITE_RESIDUAL
        //std::stringstream filename;
        //filename << "spatial_norm_"  <<"_"<< u_->time_<<"_" << u_->t_index_ << "__" << norm_counter;
        //pio_grid_function_.write(tempobject_output,filename.str().c_str());
        //auto * out = dynamic_cast<VTK_Observer<TDomain,TAlgebra>*>(out.get());
        //out->set_filename(filename.str().c_str());
        //out->step_process(tempobject_output,100000*u_->level_index+norm_counter ,u_->time,0);
        if (this->xb_out_ != SPNULL) {
            this->xb_out_->step_process(tempobject, u_->t_index_ , u_->time_,0.0, this->iteration_, 0);
        }
        //out->set_filename("output");
        norm_counter++;
        //#endif



        write_script(this->script_->SpatialNorm(u_,norm_ptr);)
        return 0;
    };*/

/*template<typename TDomain, typename TAlgebra>
int IsolatedSimpleIntegratorDriver<TDomain, TAlgebra>::Sum(double alpha, braid_Vector x_, double beta, braid_Vector y_) {
    __debug(std::cout << "GridFunctionBaseDriver::Sum" << std::endl);
    auto* xref = static_cast<SP_GridFunction *>(x_->value_);
    auto* yref = static_cast<SP_GridFunction *>(y_->value_);

    auto& xval = xref->operator*();
    auto& yval = yref->operator*();


    VecScaleAdd(yval,
                                    beta, yval,
                                    alpha, xval);

    write_script(this->script_->Sum(alpha,x_,beta,y_);)
    return 0;
};*/

/*template<typename TDomain, typename TAlgebra>
int IsolatedSimpleIntegratorDriver<TDomain, TAlgebra>::Access(braid_Vector u_, BraidAccessStatus& status) {
    __debug(std::cout << "GridFunctionBaseDriver::Access" << std::endl);
    auto ref = static_cast<SP_GridFunction *>(u_->value_)->get()->clone();

    int index;
    status.GetTIndex(&index);

    double timestamp;
    status.GetT(&timestamp);

    int iteration;
    status.GetIter(&iteration);
    this->iteration_= iteration; // todo delete

    int level;
    status.GetLevel(&level);

    int done;
    status.GetDone(&done);

    double wdt = 0;
    if (done == 1) {
        if (this->xb_out_) {
            this->xb_out_->step_process(ref, index, timestamp,wdt);
        }
        if (this->out_) {
            std::cout << " out " << std::endl;
            this->out_->step_process(ref, index, timestamp,wdt);
        }
    } else {
        if (this->xb_out_) {
            this->xb_out_->step_process(ref, index, timestamp,wdt, iteration, level);
        }
    }

    write_script(this->script_->Access(u_,status);)

    return 0;
};*/

    /*
template<typename TDomain, typename TAlgebra>
int IsolatedSimpleIntegratorDriver<TDomain, TAlgebra>::BufSize(int* size_ptr, BraidBufferStatus& status) {
    __debug(std::cout << "GridFunctionBaseDriver::BufSize" << std::endl);
    *size_ptr = 0;
#ifdef FEATURE_SPATIAL_REFINE
    *size_ptr = 0
         +sizeof(int)        // spatial-grid-level
         +sizeof(uint)       // parallel storage mask ( undefined, konsistent, unique, additive)
         +sizeof(size_t)     // number of gridfunction-elements
         +sizeof(T_VectorValueType) * (*this->u0_).size();  // size of actual vector
#else
    *size_ptr = sizeof(size_t) // number of gridfunction-elements
                 + (sizeof(T_VectorValueType) * (*this->u0).size());
    // size of actual vector

#endif


    write_script(this->script_->BufSize(size_ptr, status);)
    __debug(std::cout << "Buffer Size: " << *size_ptr << std::endl << std::flush);

    return 0;
};

template<typename TDomain, typename TAlgebra>
int IsolatedSimpleIntegratorDriver<TDomain, TAlgebra>::BufPack(braid_Vector u_, void* buffer, BraidBufferStatus& status) {
#ifdef FEATURE_SPATIAL_REFINE
    __debug(std::cout << "GridFunctionBaseDriver::BufPack" << std::endl);

    auto* u_ref = static_cast<SP_GridFunction *>(u_->value_);

    int buffer_size = 0;

    auto* chBuffer = static_cast<byte *>(buffer);
    const int spatial_level = u_ref->get()->grid_level().level();
    __debug(std::cout << "Spatial Level: " << spatial_level<< std::endl << std::flush);
    memcpy(chBuffer + buffer_size, &spatial_level, sizeof(int)); //ð
    buffer_size += sizeof(int); // ð


    uint mask = u_ref->get()->get_storage_mask(); // ð
    __debug(std::cout << "Storage Mask: " << mask<< std::endl << std::flush);
    memcpy(chBuffer + buffer_size, &mask, sizeof(uint)); //
    buffer_size += sizeof(uint); // ð

    write_script(this->script_->BufPack(u_, buffer, status,buffer_size));

    this->pack(buffer, u_ref->get(), &buffer_size);

    __debug(std::cout << "Buffer Size: " << buffer_size << std::endl << std::flush);
    __send_recv_times( std::cout << "Send t=" << timer.get() << std::endl;);
    return 0;

#else
    __debug(std::cout << "GridFunctionBaseDriver::BufPack" << std::endl);
    int buffer_size = 0; // startposition of gridfunction (will be written first) in buffer

    auto* u_ref = (SP_GridFunction*)u_->value;

    this->pack(buffer, u_ref->get(), &buffer_size);
    // buffer filled with size of vector and vector

    status.SetSize(buffer_size);

    write_script(this->script->BufPack(u_, buffer, status,buffer_size);)
    __send_recv_times( std::cout << "Send t=" << timer.get() << std::endl;);
    return 0;
#endif
};


template<typename TDomain, typename TAlgebra>
int IsolatedSimpleIntegratorDriver<TDomain, TAlgebra>::BufUnpack(void* buffer, braid_Vector* u_ptr, BraidBufferStatus& status) {
#ifdef FEATURE_SPATIAL_REFINE
            __debug(std::cout << "GridFunctionBaseDriver::BufUnpack" << std::endl);

            const auto* chBuffer = static_cast<byte *>(buffer); // ð
            int buffer_size = 0; // startposition of gridfunction (will be read first) in buffer

            auto* u = static_cast<BraidVector *>(malloc(sizeof(BraidVector)));
            *u_ptr = u;
            //ð auto* sp_u = new SP_GridFunction(new T_GridFunction(*this->u0));
            auto approx_space = this->u0_->approx_space();
            __debug(std::cout << "---------------------------- Recieved ---------------------"<< std::endl << std::flush);

            int level; // ð
            memcpy(&level, chBuffer + buffer_size, sizeof(int)); // ð
            buffer_size += sizeof(int); // ð
            __debug(std::cout << "Spatial Level: " << level<< std::endl << std::flush);

            uint mask;
            memcpy(&mask, chBuffer + buffer_size, sizeof(uint)); // ð
            buffer_size += sizeof(uint); // ð
            __debug(std::cout << "Storage Mask: " << mask<< std::endl << std::flush);



            auto* sp_u = new SP_GridFunction(new T_GridFunction(approx_space, level, false));
            sp_u->get()->set_storage_type(mask);

            write_script(this->script_->BufUnpack(buffer, u_ptr, status,buffer_size);)

            this->unpack(buffer, sp_u->get(), &buffer_size); // pos returns position of bufferpointer after writing the gridfunction
            u->value_ = sp_u;

            __debug(std::cout << "Buffer Size: " << buffer_size << std::endl << std::flush);
            __send_recv_times(std::cout << "Recv t=" << timer.get() << std::endl; );
            return 0;

#else
            __debug(std::cout << "GridFunctionBaseDriver::BufUnpack" << std::endl);
            int pos = 0; // startposition of gridfunction (will be read first) in buffer
            auto* u = (BraidVector*)malloc(sizeof(BraidVector));
            auto* sp_u = new SP_GridFunction(new T_GridFunction(*this->u0));

            this->unpack(buffer, sp_u->get(), &pos); // pos returns position of bufferpointer after writing the gridfunction
            u->value = sp_u;

            *u_ptr = u;

            write_script(this->script->BufUnpack(buffer, u_ptr, status,pos);)
            __send_recv_times(
                std::cout << "Recv t=" << timer.get() << std::endl; );
            return 0;
#endif
        };


template<typename TDomain, typename TAlgebra>
void IsolatedSimpleIntegratorDriver<TDomain, TAlgebra>::unpack(void* buffer, T_GridFunction* u_ref, int* buffer_size) {
#ifdef FEATURE_SPATIAL_REFINE
    auto* chBuffer = static_cast<byte *>(buffer);
    size_t szVector = 0;
    memcpy(&szVector, chBuffer + *buffer_size, sizeof(size_t)); // read vector size
    *buffer_size += sizeof(size_t);

    __debug(std::cout << "Recv Vector Size: " << szVector << std::endl << std::flush);

    for (size_t i = 0; i < szVector; i++) {
        T_VectorValueType val = T_VectorValueType(0);
        memcpy(&val, chBuffer + *buffer_size, sizeof(T_VectorValueType)); // read array
        *buffer_size += sizeof(T_VectorValueType);
        (*u_ref)[i] = val;
    }
    __debug(std::cout << "UnPack Buffer Position: " << buffer_size << std::endl << std::flush);
#else
    byte_t* chBuffer = (byte_t*)buffer;
    size_t szVector = 0;

    memcpy(&szVector, chBuffer + *buffer_size, sizeof(size_t)); // read vector size
    *buffer_size = sizeof(size_t);

    for (size_t i = 0; i < szVector; i++) {
        T_VectorValueType val = T_VectorValueType(0);
        memcpy(&val, chBuffer + *buffer_size, sizeof(T_VectorValueType)); // read array
        *buffer_size += sizeof(T_VectorValueType);
        (*u_ref)[i] = val;
    }
#endif


}

template<typename TDomain, typename TAlgebra>
void IsolatedSimpleIntegratorDriver<TDomain, TAlgebra>::pack(void* buffer, T_GridFunction* u_ref, int* buffer_size) {
#ifdef FEATURE_SPATIAL_REFINE
    auto* chBuffer = static_cast<byte *>(buffer);
    const size_t szVector = u_ref->size();
    __debug(std::cout << "num-elem: " << szVector << std::endl);
    memcpy(chBuffer+ *buffer_size, &szVector, sizeof(size_t)); // first value size of vector
    *buffer_size += sizeof(size_t);


    for (size_t i = 0; i < szVector; i++) {
        memcpy(chBuffer + *buffer_size, &(*u_ref)[i], sizeof(T_VectorValueType)); // array sequentially
        *buffer_size += sizeof(T_VectorValueType);
    }

#else

    byte_t* chBuffer = (byte_t*)buffer;

    size_t szVector = u_ref->size();

    memcpy(buffer, &szVector, sizeof(size_t)); // first value size of vector

    *buffer_size += sizeof(size_t);

    for (size_t i = 0; i < szVector; i++) {
        memcpy(chBuffer + *buffer_size, &(*u_ref)[i], sizeof(T_VectorValueType)); // array sequentially
        *buffer_size += sizeof(T_VectorValueType);
    }
#endif

}
*/

    /*
#ifdef FEATURE_SPATIAL_REFINE
template<typename TDomain, typename TAlgebra>
void IsolatedSimpleIntegratorDriver<TDomain, TAlgebra>:: set_level_num_ref(size_t level, int num_ref) {
        __debug(std::cout << level  << " - num ref " << num_ref << std::endl<< std::flush);
        if (this->level_num_ref.size() < level +1) {
            this->level_num_ref.resize(level + 1, 0);
        }
        this->level_num_ref[level] = num_ref;
    }
#endif* /

}}*/

#endif
