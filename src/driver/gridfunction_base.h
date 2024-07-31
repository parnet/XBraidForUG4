//todo check
#ifndef UG_PLUGIN_XBRAIDFORUG4_GRIDFUNCTION_BASE_H
#define UG_PLUGIN_XBRAIDFORUG4_GRIDFUNCTION_BASE_H

#include <ugbase.h>
#include "lib_disc/spatial_disc/domain_disc.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_algebra/vector_interface/vec_functions.h"
#include "../../../Limex/time_disc/time_integrator.hpp"
#include "../../../Limex/time_disc/linear_implicit_timestep.h"
#include "../../libs/braid/braid/braid.hpp"
#include "../core/braid_vector_struct.h"
#include "../core/space_time_communicator.h"
#include "../util/parallel_logger.h"
#include "../interface/observer_xbraid.h"
#include "../interface/initializer.h"
#include "../interface/spatial_norm.h"

namespace ug { namespace XBraidForUG4 {
    template <typename TDomain, typename TAlgebra>
    class BraidGridFunctionBase : public BraidApp {
    public:

        //--------------------------------------------------------------------------------------------------------------

        typedef typename TAlgebra::vector_type::value_type T_VectorValueType;

        typedef SmartPtr<SpaceTimeCommunicator> SP_SpaceTimeCommunicator;

        typedef IDomainDiscretization<TAlgebra> T_DomainDisc;
        typedef SmartPtr<T_DomainDisc> SP_DomainDisc;

        typedef GridFunction<TDomain, TAlgebra> T_GridFunction;
        typedef SmartPtr<T_GridFunction> SP_GridFunction;

        typedef ITimeIntegrator<TDomain, TAlgebra> T_ITimeIntegrator;
        typedef SmartPtr<T_ITimeIntegrator> SP_TimeIntegrator;

        typedef BraidInitializer<TDomain, TAlgebra> T_BraidInitializer;
        typedef SmartPtr<T_BraidInitializer> SP_BraidInitializer;

        typedef ITimeIntegratorObserver<TDomain, TAlgebra> T_ITimeIntegratorObserver;
        typedef SmartPtr<T_ITimeIntegratorObserver> SP_IObserver;

        typedef IXBraidTimeIntegratorObserver<TDomain, TAlgebra> T_IXBraidTimeIntegratorObserver;
        typedef SmartPtr<T_IXBraidTimeIntegratorObserver> SP_IXBraidTimeIntegratorObserver;

        typedef ParallelLogger T_ParallelLogger;
        typedef SmartPtr<T_ParallelLogger> SP_ParallelLogger;

        typedef BraidSpatialNorm<TDomain, TAlgebra> T_SpatialNorm;
        typedef SmartPtr<T_SpatialNorm> SP_SpatialNorm;

        //--------------------------------------------------------------------------------------------------------------

        SP_DomainDisc m_domain_disc; // used to construct ThetaTimeStep -> move?

        SP_SpaceTimeCommunicator m_comm;
        SP_ParallelLogger m_log;

        SP_GridFunction m_u0; // used for buffer and for rhs (residual) construction, todo change for adaptivity!

        SP_IObserver m_out;
        SP_IXBraidTimeIntegratorObserver m_xb_out;
        SP_BraidInitializer m_initializer;
        SP_SpatialNorm m_norm;

        bool can_residual_method = false;
        // bool can_partial_residual_method = false;
        // bool can_extrapolation = false;
        // bool can_error_estimation = false;
        // bool can_refine_timegrid = false;
        // bool can_spatial_adaptivity = false;
        // bool can_spatial_and_time_adaptivity = false;
        // bool can_tighten_conv_check = false;
        // bool can_sync = false;
        // bool can_nonlinear = false;
        // bool can_linear = false;

        bool restimate = false; // todo move? delete? change?
        bool refine_time = false;// todo move? delete? change?

        int m_levels = 0;
        int m_iteration = 0;

        //--------------------------------------------------------------------------------------------------------------

        BraidGridFunctionBase() : BraidApp(0, 0, 10, 10) {}

        BraidGridFunctionBase(MPI_Comm mpi_temporal, double tstart, double tstop, int steps) : BraidApp(mpi_temporal, tstart, tstop, steps) {}

        ~BraidGridFunctionBase() override = default;

        //--------------------------------------------------------------------------------------------------------------

        int Init(braid_Real t, braid_Vector* u_ptr) override {
            auto* u = (BraidVector*)malloc(sizeof(BraidVector));
            auto* vec = new SP_GridFunction();
            m_initializer->initialize(*vec, t);
            u->value = vec;
            u->time = t;
            *u_ptr = u;

            return 0;
        };

        int Clone(braid_Vector u_, braid_Vector* v_ptr) override {
            auto* v = (BraidVector*)malloc(sizeof(BraidVector));
            auto* uref = (SP_GridFunction*)u_->value;
            auto* vref = new SP_GridFunction();
            *vref = uref->get()->clone();
            v->value = vref;
            v->time = u_->time;
            *v_ptr = v;
            return 0;
        };

        int Free(braid_Vector u_) override {
            auto* u_value = (SP_GridFunction*)u_->value;
            delete u_value;
            free(u_);
            return 0;
        };

        int Sum(braid_Real alpha, braid_Vector x_, braid_Real beta, braid_Vector y_) override {
            auto* xref = (SP_GridFunction*)x_->value;
            auto* yref = (SP_GridFunction*)y_->value;
            ::VecAdd(beta, *yref->get(), alpha, *xref->get());
            return 0;
        };

        int SpatialNorm(braid_Vector u_, braid_Real* norm_ptr) override {
            *norm_ptr = 0;
            auto* uref = (SP_GridFunction*)u_->value;
            SP_GridFunction tempobject = uref->get()->clone(); // clone to ensure consistency
            *norm_ptr = m_norm->norm(tempobject);
            return 0;
        };

        int Access(braid_Vector u_, BraidAccessStatus& astatus) override {
            int v; // todo use v?
            int index;
            astatus.GetTIndex(&index);
            double timestamp;
            astatus.GetT(&timestamp);
            auto ref = ((SP_GridFunction*)u_->value)->get()->clone();
            int iter;
            int lvl;
            int done;
            astatus.GetIter(&iter);
            astatus.GetLevel(&lvl);
            astatus.GetDone(&done);
            double wdt = 0;
            if (done == 1) {
                if (this->m_xb_out) {
                    v = this->m_xb_out->step_process(ref, index, timestamp,wdt);
                }
                if (this->m_out) {
                    v = this->m_out->step_process(ref, index, timestamp,wdt);
                }
            } else {
                if (this->m_xb_out) {
                    v = this->m_xb_out->step_process(ref, index, timestamp,wdt, iter, lvl);
                }
            }
            return 0;
        };

        int BufSize(braid_Int* size_ptr, BraidBufferStatus& bstatus) override {
            *size_ptr = 0;
            *size_ptr = (sizeof(T_VectorValueType) * (*this->m_u0).size() // actual vector size
                    + sizeof(size_t)) // index
                + 2 * sizeof(int)
                + sizeof(double);
            return 0;
        };

        int BufPack(braid_Vector u_, void* buffer, BraidBufferStatus& bstatus) override {
            int bufferSize = 0; // startposition of gridfunction (will be written first) in buffer
            auto* u_ref = (SP_GridFunction*)u_->value;
            pack(buffer, u_ref->get(), &bufferSize); // bufferSize returns position of bufferpointer after writing the gridfunction
            char* chBuffer = (char*)buffer;
            int temprank = this->m_comm->get_temporal_rank();
            memcpy(chBuffer + bufferSize, &temprank, sizeof(int));
            bufferSize += sizeof(int);
            memcpy(chBuffer + bufferSize, &u_->time, sizeof(double));
            bufferSize += sizeof(double);
            bstatus.SetSize(bufferSize);
            return 0;
        };

        int BufUnpack(void* buffer, braid_Vector* u_ptr, BraidBufferStatus& bstatus) override {
            int pos = 0; // startposition of gridfunction (will be read first) in buffer
            auto* u = (BraidVector*)malloc(sizeof(BraidVector));
            auto* sp_u = new SP_GridFunction(new T_GridFunction(*this->m_u0));
            unpack(buffer, sp_u->get(), &pos); // pos returns position of bufferpointer after writing the gridfunction
            u->value = sp_u;
            char* chBuffer = (char*)buffer;
            int temprank;
            memcpy(&temprank, chBuffer + pos, sizeof(int));
            pos += sizeof(int);
            double t;
            memcpy(&t, chBuffer + pos, sizeof(double));
            pos += sizeof(int);
            u->time = t;
            *u_ptr = u;
            return 0;
        };

        int Sync(BraidSyncStatus& sstatus) override {
            this->m_iteration += 1;
            return 0;
        };

        int Coarsen(braid_Vector fu_, braid_Vector* cu_ptr, BraidCoarsenRefStatus& status) override {
            return 0;
        };

        int Refine(braid_Vector cu_, braid_Vector* fu_ptr, BraidCoarsenRefStatus& status) override {
            return 0;
        };

        //--------------------------------------------------------------------------------------------------------------

        void set_paralog(SP_ParallelLogger log) {
            this->m_log = log;
        }

        void init() {
            this->m_log->init();
        }

        void pack(void* buffer, T_GridFunction* u_ref, int* bufferSize) {
            char* chBuffer = (char*)buffer;
            *bufferSize = 0;
            size_t szVector = u_ref->size();
            memcpy(buffer, &szVector, sizeof(size_t));
            *bufferSize += sizeof(size_t);
            for (size_t i = 0; i < szVector; i++) {
                memcpy(chBuffer + *bufferSize, &(*u_ref)[i], sizeof(T_VectorValueType));
                *bufferSize += sizeof(T_VectorValueType);
            }
        }

        void unpack(void* buffer, T_GridFunction* u_ref, int* pos) {
            char* chBuffer = (char*)buffer;
            size_t szVector = 0;
            memcpy(&szVector, chBuffer + *pos, sizeof(size_t));
            *pos = sizeof(size_t);
            for (size_t i = 0; i < szVector; i++) {
                T_VectorValueType val = 0;
                memcpy(&val, chBuffer + *pos, sizeof(T_VectorValueType));
                *pos += sizeof(T_VectorValueType);
                (*u_ref)[i] = val;
            }
        }

        void set_time_values(double startTime, double endTime, int n) {
            this->tstart = startTime;
            this->tstop = endTime;
            this->ntime = n;
        }

        void set_start_time(double startTime) {
            this->tstart = startTime;
        }

        void set_end_time(double endTime) {
            this->tstop = endTime;
        }

        void set_number_of_timesteps(int n) {
            this->ntime = n;
        }


        void set_start_vector(SP_GridFunction p_u0) {
            this->m_u0 = p_u0;
        }

        void attach_xbraid_observer(SP_IXBraidTimeIntegratorObserver p_out) {
            this->m_xb_out = p_out;
        }

        void attach_observer(SP_IObserver p_out) {
            this->m_out = p_out;
        }

        void set_norm_provider(SP_SpatialNorm norm) {
            this->m_norm = norm;
        }

        void set_max_levels(size_t levelcount) {
            this->m_levels = levelcount;
        }

        void set_initializer(SP_BraidInitializer initializer) {
            this->m_initializer = initializer;
        }

        void set_domain(SP_DomainDisc domain) {
            this->m_domain_disc = domain;
        }

        //--------------------------------------------------------------------------------------------------------------
    };
}}
#endif //UG_PLUGIN_XBRAIDFORUG4_GRIDFUNCTION_BASE_H
