#ifndef UGPLUGIN_XBRAIDFORUG4_DRIVER_GRIDFUNCTION_BASE_HPP
#define UGPLUGIN_XBRAIDFORUG4_DRIVER_GRIDFUNCTION_BASE_HPP

#include "Limex/time_disc/time_integrator.hpp"
//2025-03 #include "../../Limex/time_disc/timestep/linear_implicit_timestep.h"
//2025-03 #include "../../Limex/time_disc/timestep/linear_implicit_timestep.h"

#include "lib_disc/function_spaces/grid_function.h"
#include "config/pragma.hpp"
#include "script_writer_driver.hpp"
#include "libs/braid/braid/braid.hpp"
#include "core/braid_vector_struct.hpp"
#include "core/space_time_communicator.hpp"
#include "util/parallel_logger.hpp"
#include "interface/observer_xbraid.hpp"
#include "interface/initializer.hpp"
#include "interface/spatial_norm.hpp"


#include "math/vector.hpp"
//2025-03  #include "../observer/vtk_observer.h"
#include <util/braid_timer.hpp>

#include "common/types.h"

#include "transfer/spatial_grid_transfer.hpp"
#include "util/parallel_io_gridfunction.hpp"




namespace ug{ namespace xbraid {

    template <typename TDomain, typename TAlgebra>
    class BraidGridFunctionBase : public BraidApp {
    public:

        //--------------------------------------------------------------------------------------------------------------

        using T_VectorValueType = typename TAlgebra::vector_type::value_type;

        using SP_SpaceTimeCommunicator = SmartPtr<SpaceTimeCommunicator> ;

        using T_DomainDisc = IDomainDiscretization<TAlgebra> ;
        using SP_DomainDisc =  SmartPtr<T_DomainDisc> ;

        using T_GridFunction=  GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;

        using T_ITimeIntegrator = ITimeIntegrator<TDomain, TAlgebra> ;
        using SP_TimeIntegrator =  SmartPtr<T_ITimeIntegrator> ;

        using T_BraidInitializer = BraidInitializer<TDomain, TAlgebra> ;
        using SP_BraidInitializer = SmartPtr<T_BraidInitializer> ;

        using T_ITimeIntegratorObserver =  ITimeIntegratorObserver<TDomain, TAlgebra> ;
        using SP_IObserver =  SmartPtr<T_ITimeIntegratorObserver> ;

        using T_IXBraidTimeIntegratorObserver =  IXBraidTimeIntegratorObserver<TDomain, TAlgebra> ;
        using SP_IXBraidTimeIntegratorObserver =  SmartPtr<T_IXBraidTimeIntegratorObserver> ;

        using T_ParallelLogger =  ParallelLogger ;
        using SP_ParallelLogger = SmartPtr<T_ParallelLogger> ;

        using T_SpatialNorm =  BraidSpatialNorm<TDomain, TAlgebra> ;
        using SP_SpatialNorm = SmartPtr<T_SpatialNorm> ;

        using T_BraidWriteScript = BraidWriteScript<TDomain, TAlgebra> ;
        using SP_BraidWriteScript = SmartPtr<T_BraidWriteScript> ;



        //--------------------------------------------------------------------------------------------------------------
    protected:
        BraidGridFunctionBase() : BraidApp(0, 0, 10, 10) {}

        BraidGridFunctionBase(MPI_Comm mpi_temporal, double tstart, double tstop, int steps) : BraidApp(mpi_temporal, tstart, tstop, steps) {}
    public:
        ~BraidGridFunctionBase() override = default;

        //--------------------------------------------------------------------------------------------------------------


        int Init(double t, braid_Vector* u_ptr) override {
            __debug(std::cout << "GridFunctionBaseDriver::Init" << std::endl);
            __send_recv_times(
                BraidTimer timer_;
                timer_.start();
                )

            auto* u = static_cast<BraidVector *>(malloc(sizeof(BraidVector)));
            auto* vec = new SP_GridFunction();
            initializer_->initialize(*vec, t);
            u->value_ = vec;
            u->time_ = t;
            *u_ptr = u;
            write_script(this->script_->Init(t, u_ptr);)

            //std::stringstream filename;
            //filename << "init_" << u->index;
            //pio_grid_function_.write(*vec,filename.str().c_str());

            /*{
                u->time = t;
                u->level_index = init_counter;
                u->level = 0;
                init_counter++ ;
            }*/
            return 0;
        };

        int Clone(braid_Vector u_, braid_Vector* v_ptr) override {
            __debug(std::cout << "GridFunctionBaseDriver::Clone" << std::endl);

            auto* v = static_cast<BraidVector *>(malloc(sizeof(BraidVector)));
            auto* uref = static_cast<SP_GridFunction *>(u_->value_);
            auto* vref = new SP_GridFunction();
            *vref = uref->get()->clone();
            v->value_ = vref;
            v->time_ = u_->time_;
            *v_ptr = v;

            /*{
                v->time = u_->time;
                v->level_index = u_->level_index;
                v->level = u_->level;
            }*/
            write_script(this->script_->Clone(u_,v_ptr);)
            return 0;
        };

        int Free(braid_Vector u_) override {
            __debug(std::cout << "GridFunctionBaseDriver::Free" << std::endl);
            write_script(this->script_->Free(u_);)

            auto* u_value = static_cast<SP_GridFunction *>(u_->value_);
            delete u_value;
            free(u_);
            return 0;
        };

// y = alpha * x + beta*y
        int Sum(double alpha, braid_Vector x_, double beta, braid_Vector y_) override {
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
        };

        int norm_counter = 0;

        int SpatialNorm(braid_Vector u_, double* norm_ptr) override {
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
        };

        int Access(braid_Vector u_, BraidAccessStatus& status) override {
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
        };

        int BufSize(int* size_ptr, BraidBufferStatus& status) override {
            __debug(std::cout << "GridFunctionBaseDriver::BufSize" << std::endl);
            *size_ptr = 0;
#ifdef FEATURE_SPATIAL_REFINE
            *size_ptr =  0
                 +sizeof(int)        // spatial-grid-level
                 +sizeof(uint)       // parallel storage mask ( undefined, konsistent, unique, additive)
                 +sizeof(size_t)     // number of gridfunction-elements
                 +sizeof(T_VectorValueType) * (*this->u0_).size();  // size of actual vector
#else
            *size_ptr =  sizeof(size_t) // number of gridfunction-elements
                         + (sizeof(T_VectorValueType) * (*this->u0).size());
            // size of actual vector

#endif


            write_script(this->script_->BufSize(size_ptr, status);)
            __debug(std::cout << "Buffer Size: " << *size_ptr << std::endl << std::flush);

            return 0;
        };

        int BufPack(braid_Vector u_, void* buffer, BraidBufferStatus& status) override {
#ifdef FEATURE_SPATIAL_REFINE
            __debug(std::cout << "GridFunctionBaseDriver::BufPack" << std::endl);

            /* unpack variables*/
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

        int BufUnpack(void* buffer, braid_Vector* u_ptr, BraidBufferStatus& status) override {
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

        int Sync(BraidSyncStatus& status) override {
            __debug(std::cout << "GridFunctionBaseDriver::Sync" << std::endl);
            this->iteration_ += 1;

            write_script(this->script_->Sync(status);)

            return 0;
        };

        void set_spatial_grid_transfer(SmartPtr<SpatialGridTransfer<TDomain,TAlgebra>> spatial_grid_transfer) {
            this->spatial_grid_transfer = spatial_grid_transfer;
        }


int Coarsen(braid_Vector           fu_,
                          braid_Vector* cu_ptr,
                          BraidCoarsenRefStatus &status) override
        {
#ifdef FEATURE_SPATIAL_REFINE
            int mgrit_level_fine;
            status.GetLevel(&mgrit_level_fine);
            // ---------------------------------------------------------- >>>


            // ---------------------------------------------------------- <<<
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
            /*{
                int t_index;
                status.GetTIndex(&t_index);
                (*cu_ptr)->time = fu_->time;
                (*cu_ptr)->level_index = t_index;
                (*cu_ptr)->level = fu_->level+1;
            }*/
            write_script(this->script_->Coarsen(fu_, cu_ptr, status);)
            __debug(std::cout << "~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ " << std::endl << std::flush);
            return 0;
#else
            this->Clone(fu_, cu_ptr); // no coarsening
            write_script(this->script->Coarsen(fu_, cu_ptr, status);)
            return 0;
#endif
        }

        /// @see braid_PtFcnSRefine.
int Refine(braid_Vector           cu_,
                                 braid_Vector          *fu_ptr,
                                 BraidCoarsenRefStatus &status) override
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
            /*{
                int t_index;
                status.GetTIndex(&t_index);
                (*fu_ptr)->time = cu_->time;
                (*fu_ptr)->level_index = t_index;
                (*fu_ptr)->level = cu_->level+1;
            }*/
            write_script(this->script_->Refine(cu_, fu_ptr, status);)
            return 0;
#else
            this->Clone(cu_, fu_ptr); // no refinement
            //void (* function)(T_GridFunction &, const T_GridFunction &));
            write_script(this->script->Coarsen(cu_, fu_ptr, status);)
            return 0;
#endif

        }



        //--------------------------------------------------------------------------------------------------------------

        void set_paralog(SP_ParallelLogger log) {
            this->log_ = log;
        }


        /**
         * after the constructor was called the init method can finish the class so that this is ready to be used
         */
        void init() {
            __send_recv_times(this->timer = BraidTimer(););
            this->log_->init();
            write_script(this->script_ = make_sp(new T_BraidWriteScript(this->comm_));) // å
        }

        void pack(void* buffer, T_GridFunction* u_ref, int* buffer_size) {
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

        void unpack(void* buffer, T_GridFunction* u_ref, int* buffer_size) {
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

        /**
         * sets initial values for the time interval and integration process like number of timesteps
         * @param startTime t_0 first timepoint in timeinterval
         * @param endTime t_N end of timeinterval
         * @param n number of timesteps
         */
        void set_time_values(double startTime, double endTime, int n) {
            this->tstart = startTime;
            this->tstop = endTime;
            this->ntime = n;
        }

        /**
         * sets the timepoint for the start of the time interval (default t_0 = 0)
         * @param startTime t_0
         */
        void set_start_time(double startTime) {
            this->tstart = startTime;
        }

        /**
         * sets the end of the timeinterval which should be integrated
         * @param endTime t_N
         */
        void set_end_time(double endTime) {
            this->tstop = endTime;
        }

        /**
         * sets the initial number of timesteps
         * @param n number of timesteps without t_0
         */
        void set_number_of_timesteps(int n) {
            this->ntime = n;
        }


        /**
         * sets a default gridfunction
         * @param p_u0 sets a u_0 for t_0
         */
        void set_start_vector(SP_GridFunction p_u0) {
            this->u0_ = p_u0;
        }

        /**
         * sets a observer which can process the intermediate results and the final results
         * @param p_out a process observer
         */
        void attach_xbraid_observer(SP_IXBraidTimeIntegratorObserver p_out) {
            this->xb_out_ = p_out;
        }

        /**
         * sets a observer which can process the results when the algorithm has finished
         * @param p_out observer object
         */
        void attach_observer(SP_IObserver p_out) {
            this->out_ = p_out;
        }

        /**
         * sets the functionality to calculate the spatial norm
         * @param norm a class method that can calculate a norm for a given gridfunction
         */
        void set_norm_provider(SP_SpatialNorm norm) {
            this->norm_ = norm;
        }

        /**
         * sets maximum number of levels
         * @param levelcount maximum number of level
         */
        void set_max_levels(size_t levelcount) {
            this->levels_ = levelcount;
        }

        /**
         * sets a generator which interpolates or clone start values for each requested timepoint
         * @param initializer initializer class
         */
        void set_initializer(SP_BraidInitializer initializer) {
            std::cout << "set-initializer" << std::endl;
            this->initializer_ = initializer;
        }

        /**
         * sets the domain for the problem
         * @param domain smart pointer to domain
         */
        void set_domain(SP_DomainDisc domain) {
            this->domain_disc_ = domain;
        }

#ifdef FEATURE_SPATIAL_REFINE
        void set_level_num_ref(size_t level, int num_ref) {
            __debug(std::cout << level  << " - num ref " << num_ref << std::endl<< std::flush);
            if (this->level_num_ref.size() < level +1) {
                this->level_num_ref.resize(level + 1, 0);
            }
            this->level_num_ref[level] = num_ref;
        }


    protected:
        std::vector<int> level_num_ref;
        SmartPtr<ApproximationSpace<TDomain>> spApproxSpace = SPNULL;
        SmartPtr<SpatialGridTransfer<TDomain,TAlgebra>> spatial_grid_transfer;

        //--------------------------------------------------------------------------------------------------------------
#else
#endif

        //--------------------------------------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------------------------------------
    public:

        //size_t m_init_counter = 0;

        SP_DomainDisc domain_disc_;
        SP_SpaceTimeCommunicator comm_;
        SP_ParallelLogger log_;

        SP_GridFunction u0_; // used for buffer and for rhs (residual) construction,

        SP_IObserver out_;
        SP_IXBraidTimeIntegratorObserver xb_out_;
        SP_BraidInitializer initializer_;
        SP_SpatialNorm norm_;
        bool can_residual_method_ = false;// error estimation and refine
        bool refine_time_ = false;// error estimation and refine
        bool restimate_ = false;
        int levels_ = 0;
        int iteration_ = 0;

        PIOGridFunction<TDomain,TAlgebra> pio_grid_function_;


        __send_recv_times(BraidTimer timer_;)
        write_script(SP_BraidWriteScript script_;)
    };

}}
#endif
