#ifndef UGPLUGIN_XBRAIDFORUG4_DRIVER_SCRIPT_WRITER_DRIVER_HPP
#define UGPLUGIN_XBRAIDFORUG4_DRIVER_SCRIPT_WRITER_DRIVER_HPP

#include "core/braid_vector_struct.hpp"
#include "util/parallel_logger.hpp"
//2025-03 #include "../core/spatial_accuracy.h"

size_t indexpool = 0;

namespace ug{ namespace xbraid {

        template<typename TDomain, typename TAlgebra>
        class BraidWriteScript  {
        public:
            using T_ParallelLogger = ParallelLogger ;
            using SP_ParallelLogger= SmartPtr<T_ParallelLogger> ;

            using T_SpaceTimeCommunicator = SpaceTimeCommunicator ;
            using SP_SpaceTimeCommunicator = SmartPtr<T_SpaceTimeCommunicator> ;

            explicit BraidWriteScript(SP_SpaceTimeCommunicator comm){
                this->comm_ = comm;

                script_log_ = make_sp(new T_ParallelLogger());
                script_log_->set_comm(comm);
                script_log_->set_filename("script");
                script_log_->init();
                script_log_->o << "initialized" << std::endl << std::flush;

            }

            ~BraidWriteScript() {
                script_log_->release();
            };

            void Init(braid_Real t, braid_Vector *u_ptr) {
                //(*u_ptr)->time = t;
                //(*u_ptr)->index = indexpool;
                //indexpool++;

                this->script_log_->o << "u_" << (*u_ptr)->index_ << " = init( time = " << t << " )" << std::endl;

            };

            void Clone(braid_Vector u_, braid_Vector *v_ptr) {
                //(*v_ptr)->index = indexpool;
                //indexpool++;

                //(*v_ptr)->time = u_->time;
                this->script_log_->o << "u_" << (*v_ptr)->index_ << " = clone( u_" << u_->index_ << " )" << std::endl;

            };

            void Free(braid_Vector u_) {
                this->script_log_->o << "u_" << u_->index_ << " = null" << std::endl;
            };

            void Sum(braid_Real alpha, braid_Vector x_, braid_Real beta, braid_Vector y_) {
                if (alpha == 0) {
                    this->script_log_->o << "u_" << y_->index_ << " = " << beta << "* u_" << y_->index_
                            << " % Skalierung "
                            <<
                            std::endl;
                } else if (beta == 0) {
                    this->script_log_->o << "u_" << y_->index_ << " = " << alpha << "*u_" << x_->index_
                            << "  % Ersetzung "
                            <<
                            std::endl;
                } else {
                    this->script_log_->o << "u_" << y_->index_ << " = " << alpha << "* u_" << x_->index_ << "  + "
                            << beta
                            << "* u_"
                            << y_->index_ << " % Summe " <<
                            std::endl;
                }
            };

            void SpatialNorm(braid_Vector u_, braid_Real *norm_ptr) {
                this->script_log_->o << "nrm = norm( u_" << u_->index_ << " )" << std::endl;
            };

            void Access(braid_Vector u_, BraidAccessStatus &status) {
                int level;
                status.GetLevel(&level);

                int num_level;
                status.GetNLevels(&num_level);

                int iteration;
                status.GetIter(&iteration);

                int t_index;
                status.GetTIndex(&t_index);

                int done;
                status.GetDone(&done);

                int num_refine;
                status.GetNRefine(&num_refine);

                int calling_function;
                status.GetCallingFunction(&calling_function);

                int delta_rank;
                status.GetDeltaRank(&delta_rank);

                int num_timepoints;
                status.GetNTPoints(&num_timepoints);

                double current_time;
                status.GetT(&current_time);

                double residual_norm;
                status.GetResidual(&residual_norm);

                double error_estimation;
                status.GetSingleErrorEstAccess(&error_estimation);

                //void GetBasisVec(braid_Vector *v_ptr, braid_Int index)

                if (done == 1) {
                    this->script_log_->o << "access( u_" << u_->index_
                                          //<< " | current_time = "<< u_->time
                                          << " ) "<< std::endl;

                } else {
                    this->script_log_->o << "access( u_" << u_->index_
                                          //<< " | current_time = "<< u_->time
                                          << " | iteration = " << iteration
                                          << " | level = " << level
                                          << " | t_index = " << t_index
                                          << " ) "<< std::endl;

                }
                this->script_log_->o << "--- "
                                          << " | level = " << level
                                          << " | num_level = "<< num_level
                                          << " | iteration = "<< iteration
                                          << " | t_index = "<< t_index
                                          << " | done = "<< done
                                          << " | num_refine = "<< num_refine
                                          << " | calling_function = "<< calling_function
                                          << " | delta_rank = "<< delta_rank
                                          << " | num_timepoints = "<< num_timepoints
                                          << " | residual_norm = "<< residual_norm
                                          << " | error_estimation = "<< error_estimation
                                          << " ) "<< std::endl;
            };

            void BufSize(braid_Int *size_ptr, BraidBufferStatus &status) {
                int message_type;
                status.GetMessageType(&message_type);
                *size_ptr += + sizeof(int) // temporal rank             ( 2Bytes )
                             + sizeof(size_t) // index                  ( 4Bytes )
                             + sizeof(double); // timestamp of the solution ( 4Bytes )

                //void SetBasisSize( braid_Int size );

                if (message_type == 0) {
                    this->script_log_->o << "sz = Step_BufferSize( )" << std::endl;
                } else if (message_type == 1) {
                    this->script_log_->o << "sz = Balance_BufferSize( )" << std::endl;
                }
            };

            void BufPack(braid_Vector u_, void *buffer, BraidBufferStatus &status, int &bufferSize) {

                int message_type;
                status.GetMessageType(&message_type);

                //bstatus.SetSize( &size );
                //bstatus.SetBasisSize( braid_Int size );


                auto* chBuffer = static_cast<byte *>(buffer);


                int temprank = this->comm_->get_temporal_rank();
                std::cout << "temporal_rank" << std::endl << std::flush;
                memcpy(chBuffer + bufferSize, &temprank, sizeof(int)); // temporal processor
                bufferSize += sizeof(int);


                size_t index = u_->index_;
                std::cout << "index" << index << std::endl << std::flush;
                memcpy(chBuffer + bufferSize, &index , sizeof(size_t)); // gridfunction timestamp
                bufferSize += sizeof(size_t);

                std::cout << "timestamp" << std::endl << std::flush;
                //double timestamp = u_->time;
                //memcpy(chBuffer + bufferSize, &timestamp, sizeof(double)); // gridfunction timestamp
                bufferSize += sizeof(double);


                if (message_type == 0) {
                    this->script_log_->o << "step_send( u_" << u_->index_ << " )" << std::endl << std::flush;
                } else if (message_type == 1) {
                    this->script_log_->o << "balance_send( u_" << u_->index_ << " )" << std::endl;
                }
                std::cout << "Writer Buffer Position: " << bufferSize  << std::endl << std::flush;
            };

            void BufUnpack(void *buffer, braid_Vector *u_ptr, BraidBufferStatus &status, int &pos) {
                int message_type;
                status.GetMessageType(&message_type);

                //bstatus.SetSize( &size );
                //bstatus.SetBasisSize( braid_Int size );
                auto* chBuffer = static_cast<byte *>(buffer);


                int temprank;
                memcpy(&temprank, chBuffer + pos, sizeof(int));
                pos += sizeof(int);

                size_t original_index;
                memcpy(&original_index, chBuffer + pos, sizeof(size_t));
                pos += sizeof(size_t);
                //(*u_ptr)->index = indexpool;
                //indexpool++;

                double timestamp;
                memcpy(&timestamp, chBuffer + pos, sizeof(double));
                pos += sizeof(double);
                //(*u_ptr)->time = timestamp;


                if (message_type == 0) {
                    this->script_log_->o << "u_" << (*u_ptr)->index_ << " = step_received( u_" << original_index<< ", proc = "<< temprank << ", t = "<< timestamp << "  )" << std::endl;
                } else if (message_type == 1) {
                    this->script_log_->o << "u_" << (*u_ptr)->index_ << " = balance_received( u _" << original_index<< ", proc = "<< temprank << ", t = "<< timestamp<<" )" << std::endl;
                }
                std::cout << "Unpack Buffer Position: " << pos  << std::endl << std::flush;
            };

            void Sync(BraidSyncStatus &status) {
                int level;
                status.GetLevel(&level);

                int num_level;
                status.GetNLevels(&num_level);

                int iteration;
                status.GetIter(&iteration);

                int num_refine;
                status.GetNRefine(&num_refine);

                int num_timepoints;
                status.GetNTPoints(&num_timepoints);

                int num_error_estimate;
                status.GetNumErrorEst(&num_error_estimate);

                int calling_function;
                status.GetCallingFunction(&calling_function);

                int done;
                status.GetDone(&done);



                this->script_log_->o  << "sync()" << std::endl;
                this->script_log_->o
                        << "---"
                        << " | level = " << level
                        << " | num_level = " << num_level
                        << " | iteration = " << iteration
                        << " | num_refine = " << num_refine
                        << " | num_timepoints = " << num_timepoints
                        << " | num_error_estimate = " << num_error_estimate
                        << " | calling_function = " << calling_function
                        << " | done = " << done
                        << std::endl << std::flush;

                double * error_estimate = (double*) malloc(sizeof(double)*num_error_estimate);
                status.GetAllErrorEst(error_estimate);
                this->script_log_->o << "=== Error Estimate" <<std::endl;
                for (int n = 0; n < num_error_estimate; n++) {
                    this->script_log_->o << "\t" << n <<"\t"<< error_estimate[n] << std::endl;
                }

                this->script_log_->o << "=== Level Configuration" <<std::endl;
                for (int l = 0; l < num_level; l++) {
                    int i_lower;
                    int i_upper;
                    status.GetTIUL(&i_upper,&i_lower,l); // input;
                    this->script_log_->o << "\t" << l << "\t" << i_lower << "\t" << i_upper << "\t timegrid = [";

                    int n_timeslice_points = i_upper - i_lower;
                    double * time_grid = (double*)malloc(sizeof(double)*n_timeslice_points);
                    status.GetTimeValues(&time_grid, i_upper, i_lower,l);
                    free(time_grid);
                    for (int i = 0; i < n_timeslice_points -1 ; i++) {
                        this->script_log_->o << time_grid[i] << " ; ";
                    }
                    this->script_log_->o << time_grid[n_timeslice_points -1];

                    this->script_log_->o << " ] " <<std::endl;

                }

                //int proc;
                //status.GetProc(&proc,in___level,in___index);

            };

            void Coarsen(braid_Vector fu_, braid_Vector *cu_ptr, BraidCoarsenRefStatus &status) {

                //(*cu_ptr)->index = indexpool;
                //indexpool++;
                //(*cu_ptr)->time = indexpool;

                int level;
                status.GetLevel(&level);

                int num_level;
                status.GetNLevels(&num_level);

                int num_refine;
                status.GetNRefine(&num_refine);

                int num_timepoints;
                status.GetNTPoints(&num_timepoints);

                int t_index;
                status.GetTIndex(&t_index);

                int iteration;
                status.GetIter(&iteration);

                double current_time;
                status.GetT(&current_time);

                double f_time_prior;
                status.GetFTprior(&f_time_prior); // fine grid time value to the left of the current time

                double f_time_stop;
                status.GetFTstop(&f_time_stop); // fine grid time value to the right// of the current time

                double c_time_prior;
                status.GetCTprior(&c_time_prior); // coarse grid time vale left of current

                double c_time_stop;
                status.GetCTstop(&c_time_stop); // coarse grid time value right of current
                ;


                this->script_log_->o << "u_" << (*cu_ptr)->index_ << " = coarsen( "
                        << " level=" << level
                        << " | t = " << current_time
                        << " | f_time_stop = " << f_time_stop
                        << " | f_time_prior = " << f_time_prior
                        << " | fu = u_" << fu_->index_
                        << " | c_time_stop = " << c_time_stop
                        << " | c_time_prior = " << c_time_prior
                << " ) " << std::endl;
                this->script_log_->o
                        << "---"
                        << " | num_level = " << num_level
                        << " | num_refine = " << num_refine
                        << " | num_timepoints = " << num_timepoints
                        << " | t_index = " << t_index
                        << " | iteration = " << iteration
                        << " ) " << std::endl << std::flush;

            };

            void Refine(braid_Vector cu_, braid_Vector *fu_ptr, BraidCoarsenRefStatus &status) {

                //(*fu_ptr)->index = indexpool;
                //indexpool++;
                // (*cu_ptr)->time = indexpool;

                int level;
                status.GetLevel(&level);

                int num_level;
                status.GetNLevels(&num_level);

                int num_refine;
                status.GetNRefine(&num_refine);

                int num_timepoints;
                status.GetNTPoints(&num_timepoints);

                int iteration;
                status.GetIter(&iteration);

                int t_index;
                status.GetTIndex(&t_index);

                double current_time;
                status.GetT(&current_time);

                double f_time_prior;
                status.GetFTprior(&f_time_prior); // fine grid time value to the left of the current time

                double f_time_stop;
                status.GetFTstop(&f_time_stop); // fine grid time value to the right// of the current time

                double c_time_prior;
                status.GetCTprior(&c_time_prior); // coarse grid time vale left of current

                double c_time_stop;
                status.GetCTstop(&c_time_stop); // coarse grid time value right of current

                this->script_log_->o << "u_" << (*fu_ptr)->index_ << " = refine( "
                        << " level=" << level
                        << " | t = " << current_time
                        << " | f_time_stop = " << f_time_stop
                        << " | f_time_prior = " << f_time_prior
                        << " | cu = u_" << cu_->index_
                        << " | c_time_stop = " << c_time_stop
                        << " | c_time_prior = " << c_time_prior
                << " ) " << std::endl;
                this->script_log_->o
                        << "---"
                        << " | num_level = " << num_level
                        << " | num_refine = " << num_refine
                        << " | num_timepoints = " << num_timepoints
                        << " | t_index = " << t_index
                        << " | iteration = " << iteration
                        << " ) " << std::endl << std::flush;
            };

            void Step(braid_Vector u_, braid_Vector ustop_, braid_Vector fstop_, BraidStepStatus &status) {
                int level;
                status.GetLevel(&level);

                int num_level;
                status.GetNLevels(&num_level);

                int done;
                status.GetDone(&done);

                int t_index;
                status.GetTIndex(&t_index);

                int iteration;
                status.GetIter(&iteration);



                int num_timepoints;
                status.GetNTPoints(&num_timepoints);

                int num_refine;
                status.GetNRefine(&num_refine);

                int calling_function;
                status.GetCallingFunction(&calling_function);

                int delta_rank;
                status.GetDeltaRank(&delta_rank);

                double tolerance;
                status.GetTol(&tolerance);

                double t_start;
                double t_stop;
                status.GetTstartTstop(&t_start, &t_stop);

                double error_estimate;
                status.GetSingleErrorEstStep(&error_estimate);

                //int rfactor;
                // int rspace;
                // status.SetRFactor(braid_Int rfactor);
                // status.SetRSpace(braid_Int rspace);
                // status.GetOldFineTolx(&old_fine_tolx_ptr);
                // status.SetOldFineTolx(braid_Real old_fine_tolx);
                // status.SetTightFineTolx(braid_Real tight_fine_tolx);
                // status.GetSpatialAccuracy(in___loose_tol, in___tight_tol, &tol_ptr); // does not used information from status / braid
                // braid_Vector v_ptr;
                // status.GetBasisVec(&v_ptr, in___index); // get vector shell from tape


                double current_dt = t_stop - t_start;



                if (fstop_ == nullptr) {
                    this->script_log_->o << "u_" << u_->index_ << " = step( "
                            << " level=" << level
                            << " | fstop = null "
                            << " | t_stop = " << t_stop
                            << " | u_stop = u_" << ustop_->index_
                            << " | t_start = " << t_start
                            << " | u = u_" << u_->index_
                            << " | dt = " << current_dt
                            << " ) " << std::endl << std::flush;

                } else {
                    this->script_log_->o << "u_" << u_->index_ << " = correcture_step( "
                            << " level=" << level
                            << " | fstop = " << fstop_->index_
                            << " | t_stop = " << t_stop
                            << " | u_stop = u_" << ustop_->index_
                            << " | t_start = " << t_start
                            << " | u = u_" << u_->index_
                            << " | dt = " << current_dt
                            << " ) " << std::endl << std::flush;
                }

                this->script_log_->o
                        << "---"
                        << " | done = " << done
                        << " | num_level = " << num_level
                        << " | num_refine = " << num_refine
                        << " | num_timepoints = " << num_timepoints
                        << " | error_estimate = " << error_estimate
                        << " | calling_function = " << calling_function
                        << " | delta_rank = " << delta_rank
                        << " | t_index = " << t_index
                        << " | iteration = " << iteration
                        << " | tolerance = " << tolerance
                        << " ---" << std::endl << std::flush;

                this->script_log_->o << "\t level | lower | upper" << std::endl;
                for (int l = 0 ; l < num_level; l++) {
                    int in_level = l;
                    int i_lower;
                    int i_upper;
                    status.GetTIUL(& i_upper, & i_lower, in_level);
                    this->script_log_->o << "\t" << l << "\t" << i_lower << "\t" << i_upper << std::endl;
                }
                this->script_log_->o << std::endl;

                int nrequest_ptr = iteration;
                double * rnorms = (double *) malloc(sizeof(double)*iteration);
                status.GetRNorms(& nrequest_ptr, rnorms);
                if (nrequest_ptr > 0) {
                    this->script_log_->o << "\t index | residual" << std::endl;
                }

                for (int n = 0 ; n < nrequest_ptr; n++) {
                    this->script_log_->o << "\t" << n << "\t" << rnorms[n] << std::endl;
                }

                if (nrequest_ptr > 0) {
                    this->script_log_->o << std::endl;
                }

                free(rnorms);


                //u_->time = t_stop;
            };

            void Residual(braid_Vector u_, braid_Vector r_, BraidStepStatus &status) {
                int level;
                status.GetLevel(&level);

                int num_level;
                status.GetNLevels(&num_level);

                int done;
                status.GetDone(&done);

                int t_index;
                status.GetTIndex(&t_index);

                int iteration;
                status.GetIter(&iteration);

                int num_timepoints;
                status.GetNTPoints(&num_timepoints);

                int num_refine;
                status.GetNRefine(&num_refine);

                int calling_function;
                status.GetCallingFunction(&calling_function);

                int delta_rank;
                status.GetDeltaRank(&delta_rank);

                double tolerance;
                status.GetTol(&tolerance);

                double t_start;
                double t_stop;
                status.GetTstartTstop(&t_start, &t_stop);

                double error_estimate;
                status.GetSingleErrorEstStep(&error_estimate);


                // int rfactor
                // status.SetRFactor(rfactor);

                // int rspace
                // status.SetRSpace(rspace);

                // double old_fine_tolx_ptr;
                // status.GetOldFineTolx(&old_fine_tolx_ptr);

                // double old_fine_tolx
                // status.SetOldFineTolx(old_fine_tolx);

                // double tight_fine_tolx
                // status.SetTightFineTolx(tight_fine_tolx);

                // double in___loose_tol;
                // double in___tight_tol;
                // double tol_ptr;
                // status.GetSpatialAccuracy(in___loose_tol, in___tight_tol, &tol_ptr); // does not used information from status / braid

                // int in___index;
                // braid_Vector v_ptr;
                // status.GetBasisVec(&v_ptr, in___index); // get vector shell from tape

                double current_dt = t_stop - t_start;


                this->script_log_->o << "u_" << r_->index_ << " = residual( "
                        << " level=" << level
                        << " | t_stop = " << t_stop
                        << " | u_stop = u_" << u_->index_
                        << " | t_start = " << t_start
                        << " | u = u_" << r_->index_
                        << " | dt = " << current_dt
                        << ") ; " << std::endl;

                this->script_log_->o
                        << "---"
                        << " | done = " << done
                        << " | num_level = " << num_level
                        << " | num_refine = " << num_refine
                        << " | num_timepoints = " << num_timepoints
                        << " | error_estimate = " << error_estimate
                        << " | calling_function = " << calling_function
                        << " | delta_rank = " << delta_rank
                        << " | t_index = " << t_index
                        << " | iteration = " << iteration
                        << " | tolerance = " << tolerance
                        << " ) " << std::endl << std::flush;

                this->script_log_->o << "\t level | lower | upper" << std::endl;
                for (int l = 0 ; l < num_level; l++) {
                    int i_upper;
                    int i_lower;
                    status.GetTIUL(& i_upper, & i_lower, l);
                    this->script_log_->o << "\t" << l << "\t" << i_lower << "\t" << i_upper << std::endl;
                }
                this->script_log_->o << std::endl;

                int nrequest_ptr;
                double * rnorms = (double *) malloc(sizeof(double)*iteration);
                nrequest_ptr = iteration;
                status.GetRNorms(& nrequest_ptr, rnorms);
                if (nrequest_ptr > 0) {
                    this->script_log_->o << "\t index | residual" << std::endl;
                }
                for (int n = 0 ; n < nrequest_ptr; n++) {
                    this->script_log_->o << "\t" << n << "\t" << rnorms[n] << std::endl;
                }
                if (nrequest_ptr > 0) {
                    this->script_log_->o << std::endl;
                }
                free(rnorms);

            }

            void set_script_log(SP_ParallelLogger log) {
                this->script_log_ = log;
            }

            SP_ParallelLogger script_log_;

            SP_SpaceTimeCommunicator comm_;


            bool write_vtk_ = false;
            bool write_gf_ = false;

            bool write_init_ = false;
            bool write_clone_ = false;
            bool write_free_ = false;

            bool write_sum_scaling_ = false;
            bool write_sum_replace_ = false;
            bool write_sum_full_ = false;
            bool write_norm_ = false;

            bool write_buf_recv_ = false;
            bool write_buf_send_ = false;

            bool write_coarsen_ = false;
            bool write_refine_ = false;



        };

    }
}
#endif