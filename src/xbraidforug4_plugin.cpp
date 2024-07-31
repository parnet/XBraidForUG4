#include "bridge/util.h"
#include "bridge/util_domain_algebra_dependent.h"

#include "core/braid_executor.h"
#include "core/space_time_communicator.h"
#include "core/replace_standard_stream.h"

#include "driver/gridfunction_base.h"
#include "driver/braid_theta_integrator_nonlinear.h"
#include "driver/braid_integrator.h"
#include "driver/braid_integrator_factory.h"
#include "driver/basic_driver.h"
#include "driver/braid_residual_stepper.h"

#include "factory/theta_integrator_factory.h"
#include "factory/fixed_step_theta_integrator_factory.h"
#include "factory/bdf_integrator_factory.h"

#include "initializer/start_value_initializer.h"
#include "initializer/zero_Initializer.h"

#include "integrator/bdf_integrator.h"
#include "integrator/theta_conststep_integrator.h"
#include "integrator/theta_conststep_integrator_nl.h"
#include "integrator/theta_integrator_nl.h"
#include "integrator/bdf_integrator_nl.h"
#include "integrator/theta_integrator.h"
#include "integrator/theta_single_timestep.h"

#include "interface/initializer.h"
#include "interface/residual_timeintegrator.h"
#include "interface/integrator_factory.h"
#include "interface/spatial_norm.h"

#include "factory/limex_integrator_factory.h"
#include "factory/simple_integrator_factory.h"
#include "factory/const_step_linear_time_integrator_factory.h"
#include "factory/linear_time_integrator_factory.h"

#include "spatial_norm/euclidian_norm.h"

#include "util/braid_timer.h"
#include "util/io_gridfunction.h"
#include "util/memory_observer.h"
#include "util/parallel_io_gridfunction.h"
#include "util/parallel_logger.h"
#include "util/world_memory.h"

#include "observer/collector_observer.h"
#include "observer/io_observer.h"
#include "observer/io_process_observer.h"
#include "observer/matlab_observer.h"
#include "observer/mock_observer.h"
#include "observer/userdata_process_evaluate_observer.h"
#include "observer/vtk_observer.h"
#include "observer/vtk_process_observer.h"
#include "observer/xb_collector_observer.h"


using namespace std;
using namespace ug::bridge;

namespace ug { namespace XBraidForUG4 {
    struct Functionality {


        template<typename TDomain, typename TAlgebra>
        static void DomainAlgebra(Registry &reg, string grp) {

            string suffix = GetDomainAlgebraSuffix<TDomain, TAlgebra>();
            string tag = GetDomainAlgebraTag<TDomain, TAlgebra>();
            /*Integrator Factories*/
            {

                {// Integrator Factory
                    typedef LimexFactory<TDomain, TAlgebra> T_LimexFactory;
                    typedef IntegratorFactory<TDomain, TAlgebra> T_IntegratorFactory;

                    string name = string("LimexFactory").append(suffix);
                    reg.add_class_<T_LimexFactory, T_IntegratorFactory>(name,grp)
                            .add_constructor()
                            .add_method("set_domain_disc", &T_LimexFactory::set_domain_disc, "", "", "")
                            .add_method("set_solver", &T_LimexFactory::set_solver, "", "", "")
                            .add_method("set_dt_max", &T_LimexFactory::set_dt_max, "", "", "")
                            .add_method("set_dt_min", &T_LimexFactory::set_dt_min, "", "", "")
                            .add_method("set_tol", &T_LimexFactory::set_tol, "", "", "")
                            .add_method("set_level_factor", &T_LimexFactory::set_level_factor, "", "", "")
                            .add_method("set_error_estimator", &T_LimexFactory::set_error_estimator, "","", "")
                            .add_method("create_time_integrator", &T_LimexFactory::create_time_integrator, "", "", "")
                            .add_method("create_leve_time_integrator", &T_LimexFactory::create_level_time_integrator,"", "", "")
                            .set_construct_as_smart_pointer(true);
                    reg.add_class_to_group(name, "LimexFactory", tag);
                }


                {
                    typedef LinearTimeIntegratorFactory<TDomain, TAlgebra> T_LinearTimeIntegratorFactory;
                    typedef IntegratorFactory<TDomain, TAlgebra> T_IntegratorFactory;
                    string name = string("LinearTimeIntegratorFactory").append(suffix);
                    reg.add_class_<T_LinearTimeIntegratorFactory, T_IntegratorFactory>(name,grp)
                            .add_constructor()
                            .add_method("set_time_disc", &T_LinearTimeIntegratorFactory::set_time_disc, "", "", "")
                            .add_method("set_solver", &T_LinearTimeIntegratorFactory::set_solver, "", "", "")
                            .add_method("create_time_integrator", &T_LinearTimeIntegratorFactory::create_time_integrator,"", "", "")
                            .add_method("create_level_time_integrator",&T_LinearTimeIntegratorFactory::create_level_time_integrator, "", "", "")
                            .set_construct_as_smart_pointer(true);
                    reg.add_class_to_group(name, "LinearTimeIntegratorFactory", tag);
                }

                {
                    typedef ConstStepLinearTimeIntegratorFactory<TDomain, TAlgebra> T_ConstStepLinearTimeIntegratorFactory;
                    typedef IntegratorFactory<TDomain, TAlgebra> T_IntegratorFactory;
                    string name = string("ConstStepLinearTimeIntegratorFactory").append(suffix);
                    reg.add_class_<T_ConstStepLinearTimeIntegratorFactory, T_IntegratorFactory>(name,grp)
                            .add_constructor()
                            .add_method("set_time_disc", &T_ConstStepLinearTimeIntegratorFactory::set_time_disc, "", "", "")
                            .add_method("set_solver", &T_ConstStepLinearTimeIntegratorFactory::set_solver, "", "", "")
                            .add_method("create_time_integrator", &T_ConstStepLinearTimeIntegratorFactory::create_time_integrator,"", "", "")
                            .add_method("create_level_time_integrator",&T_ConstStepLinearTimeIntegratorFactory::create_level_time_integrator, "", "", "")
                            .add_method("set_num_steps", &T_ConstStepLinearTimeIntegratorFactory::set_num_steps, "", "", "")
                            .set_construct_as_smart_pointer(true);
                    reg.add_class_to_group(name, "ConstStepLinearTimeIntegratorFactory", tag);
                }
            }

            /*Integrator Methods*/
            {
                {
                    typedef BDF_Integrator<TDomain, TAlgebra> T_BDF_Integrator;
                    typedef ITimeIntegrator<TDomain, TAlgebra> TBase;
                    string name = string("BDF_Integrator").append(suffix);
                    reg.add_class_<T_BDF_Integrator, TBase>(name, grp)
                            .add_constructor()
                            .add_method("set_domain", &T_BDF_Integrator::set_domain, "", "", "")
                            .add_method("set_solver", &T_BDF_Integrator::set_solver, "", "", "")
                            .add_method("set_order", &T_BDF_Integrator::set_order, "", "", "")
                            .add_method("set_reassemble_threshold", &T_BDF_Integrator::set_reassemble_threshold, "","", "")
                            .set_construct_as_smart_pointer(true);
                    reg.add_class_to_group(name, "BDF_Integrator", tag);
                }
                {
                    typedef BDF_IntegratorNL<TDomain, TAlgebra> T_NLBDFIntegrator;
                    typedef ITimeIntegrator<TDomain, TAlgebra> T_ITimeIntegrator;
                    string name = string("BDF_IntegratorNL").append(suffix);
                    reg.add_class_<T_NLBDFIntegrator, T_ITimeIntegrator>(name, grp)
                            .add_constructor()
                            .add_method("set_domain", &T_NLBDFIntegrator::set_domain, "", "", "")
                            .add_method("set_solver", &T_NLBDFIntegrator::set_solver, "", "", "")
                            .add_method("set_order", &T_NLBDFIntegrator::set_order, "", "", "")
                            .set_construct_as_smart_pointer(true);
                    reg.add_class_to_group(name, "BDF_IntegratorNL", tag);
                }
                 {
                typedef ThetaIntegrator<TDomain, TAlgebra> ThetaIntegrator;
                typedef ITimeIntegrator<TDomain, TAlgebra> T_Base;
                string name = string("ThetaIntegrator").append(suffix);
                reg.add_class_<ThetaIntegrator, T_Base>(name, grp)
                        .add_constructor()
                        .add_method("set_domain", &ThetaIntegrator::set_domain, "", "", "")
                        .add_method("set_solver", &ThetaIntegrator::set_solver, "", "", "")
                        .add_method("set_theta", &ThetaIntegrator::set_theta, "", "", "")
                        .add_method("set_reassemble_threshold", &ThetaIntegrator::set_reassemble_threshold, "","", "")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name, "ThetaIntegrator", tag);
            }                 {
                typedef ThetaConstStepIntegrator<TDomain, TAlgebra> T_ThetaConstStepIntegrator;
                typedef ITimeIntegrator<TDomain, TAlgebra> T_Base;
                string name = string("ThetaConstStepIntegrator").append(suffix);
                reg.add_class_<T_ThetaConstStepIntegrator, T_Base>(name, grp)
                        .add_constructor()
                        .add_method("set_domain", &T_ThetaConstStepIntegrator::set_domain, "", "", "")
                        .add_method("set_solver", &T_ThetaConstStepIntegrator::set_solver, "", "", "")
                        .add_method("set_theta", &T_ThetaConstStepIntegrator::set_theta, "", "", "")
                        .add_method("set_num_steps", &T_ThetaConstStepIntegrator::set_num_steps, "", "", "")
                        .add_method("set_reassemble_threshold", &T_ThetaConstStepIntegrator::set_reassemble_threshold, "","", "")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name, "ThetaConstStepIntegrator", tag);
            }
            {
                typedef ThetaIntegratorNL<TDomain, TAlgebra> T_ThetaIntegratorNL;
                typedef ITimeIntegrator<TDomain, TAlgebra> T_ITimeIntegrator;
                string name = string("ThetaIntegratorNL").append(suffix);
                reg.add_class_<T_ThetaIntegratorNL, T_ITimeIntegrator>(name, grp)
                        .add_constructor()
                        .add_method("set_domain", &T_ThetaIntegratorNL::set_domain, "", "", "")
                        .add_method("set_solver", &T_ThetaIntegratorNL::set_solver, "", "", "")
                        .add_method("set_theta", &T_ThetaIntegratorNL::set_theta, "", "", "")
                        .add_method("apply", &T_ThetaIntegratorNL::apply, "", "", "")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name, "ThetaIntegratorNL", tag);
            }
            {
                typedef ThetaConstStepIntegratorNL<TDomain, TAlgebra> T_ThetaConstStepIntegratorNL;
                typedef ITimeIntegrator<TDomain, TAlgebra> T_ITimeIntegrator;
                string name = string("ThetaConstStepIntegratorNL").append(suffix);
                reg.add_class_<T_ThetaConstStepIntegratorNL, T_ITimeIntegrator>(name, grp)
                        .add_constructor()
                        .add_method("set_domain", &T_ThetaConstStepIntegratorNL::set_domain, "", "", "")
                        .add_method("set_solver", &T_ThetaConstStepIntegratorNL::set_solver, "", "", "")
                        .add_method("set_theta", &T_ThetaConstStepIntegratorNL::set_theta, "", "", "")
                        .add_method("set_num_steps", &T_ThetaConstStepIntegratorNL::set_num_steps, "", "", "")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name, "ThetaConstStepIntegratorNL", tag);
            }



            {
                typedef IResidualTimeIntegrator<TDomain, TAlgebra> T_IResidualTimeIntegrator;
                string name = string("IResidualTimeIntegrator").append(suffix);
                reg.add_class_<T_IResidualTimeIntegrator>(name, grp);
                reg.add_class_to_group(name, "IResidualTimeIntegrator", tag);
            }
            {
                typedef ThetaSingleTimeStep<TDomain, TAlgebra> T_ThetaSingleTimeStep;
                typedef IResidualTimeIntegrator<TDomain, TAlgebra> T_IResidualTimeIntegrator;
                string name = string("ThetaSingleTimeStep").append(suffix);
                reg.add_class_<T_ThetaSingleTimeStep, T_IResidualTimeIntegrator>(name, grp)
                        .add_constructor()
                        .add_method("set_domain", &T_ThetaSingleTimeStep::set_domain, "", "", "")
                        .add_method("set_solver", &T_ThetaSingleTimeStep::set_solver, "", "", "")
                        .add_method("set_theta", &T_ThetaSingleTimeStep::set_theta, "", "", "")
                        .add_method("set_reassemble_threshold", &T_ThetaSingleTimeStep::set_reassemble_threshold, "","", "")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name, "ThetaSingleTimeStep", tag);
            }
            }

            /* Utility */
            {
                {
                    typedef IOGridFunction<TDomain, TAlgebra> T_IOGridFunction;
                    string name = string("IOGridFunction").append(suffix);
                    reg.add_class_<T_IOGridFunction>(name, grp)
                            .add_constructor()
                            .add_method("write", &T_IOGridFunction::write, "", "", "")
                            .add_method("read", &T_IOGridFunction::read, "", "", "")
                            .set_construct_as_smart_pointer(true);
                    reg.add_class_to_group(name, "IOGridFunction", tag);
                }
                {
                    typedef PIOGridFunction<TDomain, TAlgebra> T_PIOGridFunction;
                    string name = string("PIOGridFunction").append(suffix);
                    reg.add_class_<T_PIOGridFunction>(name, grp)
                            .add_constructor()
                            .add_method("write", &T_PIOGridFunction::write, "", "", "")
                            .add_method("read", &T_PIOGridFunction::read, "", "", "")
                            .set_construct_as_smart_pointer(true);
                    reg.add_class_to_group(name, "PIOGridFunction", tag);
                }
            }




            {// BraidInitializer
                typedef BraidInitializer<TDomain, TAlgebra> T_BraidInitializer;
                typedef GridFunction<TDomain, TAlgebra> T_GridFunction;
                typedef SmartPtr<T_GridFunction> SP_GridFunction;

                string name = string("BraidInitializer").append(suffix);
                reg.add_class_<T_BraidInitializer>(name, grp)
                    //.add_method("initialize", static_cast<void (T_BraidInitializer::*)(SP_GridFunction & ,number)>(&T_BraidInitializer::initialize), "", "", "")
                    .add_method("set_start_values", &T_BraidInitializer::set_start_values, "", "", "")
                    .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name, "BraidInitializer", tag);
            }

            {//start_value_initializer
                typedef GridFunctionInitializer<TDomain, TAlgebra> T_StartValueInitializer;
                typedef BraidInitializer<TDomain, TAlgebra> T_BraidInitializer;
                string name = string("GridFunctionInitializer").append(suffix);
                reg.add_class_<T_StartValueInitializer, T_BraidInitializer>(name, grp)
                        .add_constructor()
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name, "GridFunctionInitializer", tag);
            }

            {//ZeroInitializer
                typedef ZeroInitializer<TDomain, TAlgebra> T_ZeroValueInitializer;
                typedef BraidInitializer<TDomain, TAlgebra> T_BraidInitializer;
                string name = string("ZeroInitializer").append(suffix);
                reg.add_class_<T_ZeroValueInitializer, T_BraidInitializer>(name, grp)
                        .add_constructor()
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name, "ZeroInitializer", tag);
            }

            { // Observer Interface
                typedef IXBraidTimeIntegratorObserver<TDomain, TAlgebra> T_IXBraidTimeIntegratorObserver;
                typedef ITimeIntegratorObserver<TDomain, TAlgebra> T_ITimeIntegratorObserver;
                string name = string("XBraidTimeIntegratorObserver").append(suffix);
                reg.add_class_<T_IXBraidTimeIntegratorObserver, T_ITimeIntegratorObserver>(name, grp)
                        .add_method("write", &T_IXBraidTimeIntegratorObserver::write, "", "","")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name, "XBraidTimeIntegratorObserver", tag);
            }

            {// Normal Collector Depcrated by Integrator Subject
                typedef TimeIntegratorObserverCollector<TDomain, TAlgebra> T_TimeIntegratorObserverCollector;
                typedef ITimeIntegratorObserver<TDomain, TAlgebra> T_ITimeIntegratorObserver;
                string name_multi = string("TimeIntegratorObserverCollector").append(suffix);
                reg.add_class_<T_TimeIntegratorObserverCollector, T_ITimeIntegratorObserver>(name_multi, grp)
                        .add_constructor()
                        .add_method("attach_observer", &T_TimeIntegratorObserverCollector::attach_observer, "", "", "")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name_multi, "TimeIntegratorObserverCollector", tag);
            }

            {// XBraid Observer Collector
                typedef XBraidTimeIntegratorObserverCollector<TDomain, TAlgebra> T_XBraidTimeIntegratorObserverCollector;
                typedef IXBraidTimeIntegratorObserver<TDomain, TAlgebra> T_IXBraidTimeIntegratorObserver;
                typedef GridFunction<TDomain, TAlgebra> T_GridFunction;
                typedef SmartPtr<T_GridFunction> SP_GridFunction;
                string name_multi = string("XBraidTimeIntegratorObserverCollector").append(suffix);
                reg.add_class_<T_XBraidTimeIntegratorObserverCollector, T_IXBraidTimeIntegratorObserver>(name_multi, grp)
                        .add_constructor()
                        .add_method("attach_observer", &T_XBraidTimeIntegratorObserverCollector::attach_observer, "", "", "")
                        .add_method("attach_common_observer", &T_XBraidTimeIntegratorObserverCollector::attach_common_observer, "", "", "")
                        .add_method("step_process",(bool (T_XBraidTimeIntegratorObserverCollector::*)(SP_GridFunction u, int, number, number)) &T_XBraidTimeIntegratorObserverCollector::step_process, "", "", "")
                        .add_method("step_process", (bool (T_XBraidTimeIntegratorObserverCollector::*)(SP_GridFunction u, int, number, number, int, int)) &T_XBraidTimeIntegratorObserverCollector::step_process, "","", "")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name_multi, "XBraidTimeIntegratorObserverCollector", tag);
            }

            {// Normal VTK Observer
                typedef VTK_Observer<TDomain, TAlgebra> T_VTK_Observer;
                typedef ITimeIntegratorObserver<TDomain, TAlgebra> T_ITimeIntegratorObserver;
                typedef VTKOutput<TDomain::dim> T_VTKOutput;
                typedef SmartPtr<T_VTKOutput> SP_VTKOutput;
                string name_multi = string("VTK_Observer").append(suffix);
                reg.add_class_<T_VTK_Observer, T_ITimeIntegratorObserver>(name_multi, grp)
                        .template add_constructor<void (*)(SP_VTKOutput, const char *)>("")
                        .add_method("step_process", &T_VTK_Observer::step_process, "", "", "")
                        .add_method("write_time_pvd", &T_VTK_Observer::write_time_pvd, "", "", "")
                        .add_method("set_filename", &T_VTK_Observer::set_filename, "", "", "")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name_multi, "TimeIntegratorObserverCollector", tag);
            }



            {//
                typedef MATLAB_Observer<TDomain, TAlgebra> T_MATLAB_Observer;
                typedef ITimeIntegratorObserver<TDomain, TAlgebra> T_ITimeIntegratorObserver;

                typedef GridFunction<TDomain, TAlgebra> T_GridFunction;
                typedef SmartPtr<T_GridFunction> SP_GridFunction;


                typedef SmartPtr<ParallelLogger> SP_ParallelLogger;
                string name_multi = string("MATLAB_Observer").append(suffix);
                reg.add_class_<T_MATLAB_Observer, T_ITimeIntegratorObserver>(name_multi, grp)
                .template add_constructor<void (*)(SP_ParallelLogger)>(", ")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name_multi, "MATLAB_Observer", tag);
            }

            {// XBraid VTK Process Observer
                typedef VTK_ProcessObserver<TDomain, TAlgebra> T_VTK_ProcessObserver;
                typedef IXBraidTimeIntegratorObserver<TDomain, TAlgebra> T_IXBraidTimeIntegratorObserver;

                typedef GridFunction<TDomain, TAlgebra> T_GridFunction;
                typedef SmartPtr<T_GridFunction> SP_GridFunction;

                typedef VTKOutput<TDomain::dim> T_VTKOutput;
                typedef SmartPtr<T_VTKOutput> SP_VTKOutput;

                string name_multi = string("VTK_ProcessObserver").append(suffix);
                reg.add_class_<T_VTK_ProcessObserver, T_IXBraidTimeIntegratorObserver>(name_multi, grp)
                        .template add_constructor<void (*)(SP_VTKOutput, const char *)>("")
                        .add_method("step_process", (bool (T_VTK_ProcessObserver::*)(SP_GridFunction u, int,number,number) ) &T_VTK_ProcessObserver::step_process, "","", "")
                        .add_method("step_process", (bool (T_VTK_ProcessObserver::*)(SP_GridFunction u, int,number,number, int,int) ) &T_VTK_ProcessObserver::step_process, "","", "")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name_multi, "VTK_ProcessObserver", tag);
            }

            {// EvalScriptor
                typedef UserdataProcessEvaluateObserver<TDomain, TAlgebra> T_UserdataProcessEvaluateObserver;
                typedef IXBraidTimeIntegratorObserver<TDomain, TAlgebra> T_IXBraidTimeIntegratorObserver;
                string name_eval = string("EvalScriptor").append(suffix);
                reg.add_class_<T_UserdataProcessEvaluateObserver, T_IXBraidTimeIntegratorObserver>(name_eval, grp)
                        .add_constructor()
                        .add_method("setFile",static_cast<void (T_UserdataProcessEvaluateObserver::*)(const char *)>(&T_UserdataProcessEvaluateObserver::set_file),"","", "")
                        .add_method("setGeneratorComponent", &T_UserdataProcessEvaluateObserver::set_generator_component, "","", "")
                        .add_method("setVectorGenerator", static_cast<void (T_UserdataProcessEvaluateObserver::*)(const char *)>(&T_UserdataProcessEvaluateObserver::set_vector_generator), "", "","")
                        .add_method("setDomain", &T_UserdataProcessEvaluateObserver::set_domain, "", "","")
                        .add_method("setRelative", &T_UserdataProcessEvaluateObserver::set_relative, "","", "")
                        .add_method("write_time_pvd", &T_UserdataProcessEvaluateObserver::write_time_pvd, "", "","")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name_eval, "EvalScriptor", tag);
            }

            {//Norm Provider
                typedef BraidSpatialNorm<TDomain, TAlgebra> T_BraidSpatialNorm;
                string name = string("BraidSpatialNorm").append(suffix);
                reg.add_class_<T_BraidSpatialNorm>(name, grp);
                reg.add_class_to_group(name, "BraidSpatialNorm", tag);
            }

            {
                typedef BraidEuclidianNorm<TDomain, TAlgebra> T_BraidEuclidianNorm;
                typedef BraidSpatialNorm<TDomain, TAlgebra> T_BraidSpatialNorm;
                string name = string("BraidEuclidianNorm").append(suffix);
                reg.add_class_<T_BraidEuclidianNorm, T_BraidSpatialNorm>(name, grp)
                        .add_constructor()
                        .add_method("norm", &T_BraidEuclidianNorm::norm, "", "", "")
                        .set_construct_as_smart_pointer(true);;
                reg.add_class_to_group(name, "BraidEuclidianNorm", tag);
            }


            {//SimpleLimexFactory
                typedef IntegratorFactory<TDomain, TAlgebra> T_IntegratorFactory;
                string name = string("IntegratorFactory").append(suffix);
                reg.add_class_<T_IntegratorFactory>(name, grp);
                reg.add_class_to_group(name, "IntegratorFactory", tag);
            }
            { // Theta Factory
                typedef ThetaIntegratorFactory<TDomain, TAlgebra> T_ThetaIntegrator;
                typedef IntegratorFactory<TDomain, TAlgebra> T_IntegratorFactory;
                string name = string("ThetaIntegratorFactory").append(suffix);
                reg.add_class_<T_ThetaIntegrator, T_IntegratorFactory>(name, grp)
                        .add_constructor()
                        .add_method("set_domain", &T_ThetaIntegrator::set_domain, "", "","")
                        .add_method("set_solver", &T_ThetaIntegrator::set_solver, "", "", "")
                        .add_method("set_theta", &T_ThetaIntegrator::set_theta, "", "", "")
                        .add_method("set_level_theta", &T_ThetaIntegrator::set_level_theta, "", "", "")
                        .add_method("create_time_integrator", &T_ThetaIntegrator::create_time_integrator, "", "","")
                        .add_method("create_level_time_integrator",&T_ThetaIntegrator::create_level_time_integrator, "", "", "")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name, "ThetaIntegratorFactory", tag);
            }
            {
                typedef FixedStepThetaIntegratorFactory<TDomain, TAlgebra> T_FixedStepThetaIntegratorFactory;
                typedef IntegratorFactory<TDomain, TAlgebra> TBase;
                string name = string("FixedStepThetaIntegratorFactory").append(suffix);
                reg.add_class_<T_FixedStepThetaIntegratorFactory, TBase>(name, grp)
                        .add_constructor()
                        .add_method("set_domain", &T_FixedStepThetaIntegratorFactory::set_domain, "", "", "")
                        .add_method("set_solver", &T_FixedStepThetaIntegratorFactory::set_solver, "", "", "")
                        .add_method("set_theta", &T_FixedStepThetaIntegratorFactory::set_theta, "", "", "")
                        .add_method("set_level_theta", &T_FixedStepThetaIntegratorFactory::set_level_theta, "", "", "")
                        .add_method("set_level_num_steps", &T_FixedStepThetaIntegratorFactory::set_level_num_steps, "", "", "")
                        .add_method("create_time_integrator", &T_FixedStepThetaIntegratorFactory::create_time_integrator, "", "", "")
                        .add_method("create_level_time_integrator", &T_FixedStepThetaIntegratorFactory::create_level_time_integrator, "", "", "")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name, "FixedStepThetaIntegratorFactory", tag);
            }
            {// Theta Factory
                typedef BDF_IntegratorFactory<TDomain, TAlgebra> T_BDF_IntegratorFactory;
                typedef IntegratorFactory<TDomain, TAlgebra> TBase;
                string name = string("BDF_IntegratorFactory").append(suffix);
                reg.add_class_<T_BDF_IntegratorFactory, TBase>(name, grp)
                        .add_constructor()
                        .add_method("set_domain", &T_BDF_IntegratorFactory::set_domain, "", "", "")
                        .add_method("set_solver", &T_BDF_IntegratorFactory::set_solver, "", "", "")
                        .add_method("set_order", &T_BDF_IntegratorFactory::set_order, "", "", "")
                        .add_method("set_level_order", &T_BDF_IntegratorFactory::set_level_order, "", "", "")
                        .add_method("create_time_integrator", &T_BDF_IntegratorFactory::create_time_integrator, "", "", "")
                        .add_method("create_level_time_integrator", &T_BDF_IntegratorFactory::create_level_time_integrator, "", "", "")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name, "BDF_IntegratorFactory", tag);
            }
            {// SimpleIntegratorFactory
                typedef SimpleIntegratorFactory<TDomain, TAlgebra> T_SimpleIntegratorFactory;
                typedef IntegratorFactory<TDomain, TAlgebra> T_IntegratorFactory;
                string name = string("SimpleIntegratorFactory").append(suffix);
                reg.add_class_<T_SimpleIntegratorFactory, T_IntegratorFactory>(name,grp)
                        .add_constructor()
                        .add_method("set_domain_disc", &T_SimpleIntegratorFactory::set_domain_disc, "", "", "")
                        .add_method("set_solver", &T_SimpleIntegratorFactory::set_solver, "", "", "")
                        .add_method("set_dt_min", &T_SimpleIntegratorFactory::set_dt_min, "", "end time", "")
                        .add_method("set_dt_max", &T_SimpleIntegratorFactory::set_dt_max, "", "", "")
                        .add_method("create_time_integrator", &T_SimpleIntegratorFactory::create_time_integrator,"", "", "")
                        .add_method("create_level_time_integrator",&T_SimpleIntegratorFactory::create_level_time_integrator, "", "", "")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name, "SimpleIntegratorFactory", tag);
            }

            {// Braid Grid Function Base
                typedef BraidGridFunctionBase<TDomain, TAlgebra> T_BraidGridFunctionBase;
                string name = string("BraidGridFunctionBase").append(suffix);
                reg.add_class_<T_BraidGridFunctionBase>(name, grp)
                        .add_method("init", &T_BraidGridFunctionBase::init, "", "", "")
                        .add_method("set_start_time", &T_BraidGridFunctionBase::set_start_time, "", "")
                        .add_method("set_end_time", &T_BraidGridFunctionBase::set_end_time, "", "", "")
                        .add_method("set_number_of_timesteps", &T_BraidGridFunctionBase::set_number_of_timesteps, "", "", "")
                        .add_method("set_time_values", static_cast<void (T_BraidGridFunctionBase::*)(double, double, int)>(&T_BraidGridFunctionBase::set_time_values), "", "", "")
                        .add_method("set_start_vector", &T_BraidGridFunctionBase::set_start_vector, "", "", "")
                        .add_method("set_norm_provider", &T_BraidGridFunctionBase::set_norm_provider, "", "", "")
                        .add_method("attach_xbraid_observer", &T_BraidGridFunctionBase::attach_xbraid_observer, "", "", "")
                        .add_method("attach_observer", &T_BraidGridFunctionBase::attach_observer, "", "", "")
                        .add_method("set_max_levels", &T_BraidGridFunctionBase::set_max_levels, "", "")
                        .add_method("set_domain", &T_BraidGridFunctionBase::set_domain, "", "");
                reg.add_class_to_group(name, "BraidGridFunctionBase", tag);
            }
            { // Braid Time Integrator
                typedef BraidIntegrator<TDomain, TAlgebra> T_BraidIntegrator;
                typedef BraidGridFunctionBase<TDomain, TAlgebra> T_BraidGridFunctionBase;
                string name = string("BraidIntegrator").append(suffix);
                reg.add_class_<T_BraidIntegrator, T_BraidGridFunctionBase>(name, grp)
                        .add_constructor()
                        .add_method("print_settings", &T_BraidIntegrator::print_settings, "", "", "")
                        .add_method("set_default_integrator", &T_BraidIntegrator::set_default_integrator, "", "", "")
                        .add_method("set_ref_factor", &T_BraidIntegrator::set_ref_factor, "", "", "")
                        .add_method("set_threshold", &T_BraidIntegrator::set_threshold, "", "", "")
                        .add_method("set_integrator", &T_BraidIntegrator::set_integrator, "", "", "")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name, "BraidIntegrator", tag);
            }
            {
                typedef BraidNLIntegrator<TDomain, TAlgebra> T_BraidNLIntegrator;
                typedef BraidGridFunctionBase<TDomain, TAlgebra> T_BraidGridFunctionBase;
                string name = string("BraidNLIntegrator").append(suffix);
                reg.add_class_<T_BraidNLIntegrator, T_BraidGridFunctionBase>(name, grp)
                        .add_constructor()
                        .add_method("print_settings", &T_BraidNLIntegrator::print_settings, "", "", "")
                        .add_method("set_default_integrator", &T_BraidNLIntegrator::set_default_integrator, "", "", "")
                        .add_method("set_ref_factor", &T_BraidNLIntegrator::set_ref_factor, "", "", "")
                        .add_method("set_threshold", &T_BraidNLIntegrator::set_threshold, "", "", "")
                        .add_method("set_integrator", &T_BraidNLIntegrator::set_integrator, "", "", "")
                        .add_method("set_conv_check", &T_BraidNLIntegrator::set_conv_check, "", "", "")
                        .add_method("set_tol", &T_BraidNLIntegrator::set_tol, "", "", "")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name, "BraidNLIntegrator", tag);
            }
            {// Braid Time Integrator Factory
                typedef BraidIntegratorFactory<TDomain, TAlgebra> T_BraidIntegratorFactory;
                typedef BraidGridFunctionBase<TDomain, TAlgebra> T_BraidGridFunctionBase;
                string name = string("BraidIntegratorFactory").append(suffix);
                reg.add_class_<T_BraidIntegratorFactory, T_BraidGridFunctionBase>(name, grp)
                        .add_constructor()
                        .add_method("print_settings", &T_BraidIntegratorFactory::print_settings, "", "", "")
                        .add_method("set_fine_time_integrator", &T_BraidIntegratorFactory::set_fine_time_integrator, "", "", "")
                        .add_method("set_coarse_time_integrator", &T_BraidIntegratorFactory::set_coarse_time_integrator, "", "", "")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name, "BraidIntegratorFactory", tag);
            }
            {// Braid Time Stepper
                typedef BasicDriver<TDomain, TAlgebra> T_BasicDriver;
                typedef BraidGridFunctionBase<TDomain, TAlgebra> T_BraidGridFunctionBase;
                string name = string("BasicDriver").append(suffix);
                reg.add_class_<T_BasicDriver, T_BraidGridFunctionBase>(name, grp)
                        .add_constructor()
                        .add_method("print_settings", &T_BasicDriver::print_settings, "", "", "")
                        .add_method("set_domain", &T_BasicDriver::set_domain, "", "", "")
                        .add_method("set_default_integrator", &T_BasicDriver::set_default_integrator, "", "","")
                        .add_method("set_integrator", &T_BasicDriver::set_integrator, "", "", "")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name, "BasicDriver", tag);
            }
            {// Braid Time Stepper
                typedef BraidResidualStepper<TDomain, TAlgebra> T_BraidTimeStepper;
                typedef BraidGridFunctionBase<TDomain, TAlgebra> T_BraidGridFunctionBase;
                string name = string("BraidResidualStepper").append(suffix);
                reg.add_class_<T_BraidTimeStepper, T_BraidGridFunctionBase>(name, grp)
                        .add_constructor()
                        .add_method("print_settings", &T_BraidTimeStepper::print_settings, "", "", "")
                        .add_method("set_approx_space", &T_BraidTimeStepper::set_approx_space, "", "", "")
                        .add_method("set_adapt_convergence", &T_BraidTimeStepper::setAdaptConv, "", "")
                        .add_method("set_domain", &T_BraidTimeStepper::set_domain, "", "")
                        .add_method("set_solver", &T_BraidTimeStepper::set_solver, "", "")
                        .add_method("set_force_convergence", &T_BraidTimeStepper::setForceConv, "", "", "")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name, "BraidResidualStepper", tag);
            }// Grid Function Executor
            {
                typedef BraidExecutor<TDomain, TAlgebra> T_BraidExecutor;
                typedef SmartPtr<SpaceTimeCommunicator> SP_SpaceTimeCommunicator;
                typedef BraidGridFunctionBase<TDomain, TAlgebra> T_BraidGridFunctionBase;
                typedef SmartPtr<T_BraidGridFunctionBase> SP_BraidGridFunctionBase;
                typedef GridFunction<TDomain, TAlgebra> T_GridFunction;
                typedef SmartPtr<T_GridFunction> SP_GridFunction;
                string name = string("BraidExecutor").append(suffix);
                reg.add_class_<T_BraidExecutor>(name, grp)
                        .template add_constructor<void (*)(SP_SpaceTimeCommunicator, SP_BraidGridFunctionBase)>(", ")
                        .add_method("apply", static_cast<bool (T_BraidExecutor::*)(SP_GridFunction, number, SP_GridFunction, number)>(&T_BraidExecutor::apply),"", "", "")
                        .add_method("set_residual", &T_BraidExecutor::set_residual, "", "", "")
                        .add_method("set_n_relax", &T_BraidExecutor::set_n_relax, "", "", "")
                        .add_method("set_c_factor", &T_BraidExecutor::set_c_factor, "", "", "")
                        .add_method("set_max_levels", &T_BraidExecutor::set_max_levels, "", "", "")
                        .add_method("set_skip_downcycle_work", &T_BraidExecutor::set_skip_downcycle_work, "", "", "")
                        .add_method("set_min_coarse", &T_BraidExecutor::set_min_coarse, "", "", "")
                        .add_method("set_max_iterations", &T_BraidExecutor::set_max_iterations, "", "", "")
                        .add_method("set_absolute_tol", &T_BraidExecutor::set_absolute_tol, "", "", "")
                        .add_method("set_relative_tol", &T_BraidExecutor::set_relative_tol, "", "", "")
                        .add_method("set_temporal_norm", &T_BraidExecutor::set_temporal_norm, "", "", "")
                        .add_method("set_sequential", &T_BraidExecutor::set_sequential, "", "", "")
                        .add_method("set_store_values", &T_BraidExecutor::set_store_values, "", "", "")
                        .add_method("set_spatial_coarsen_and_refine", &T_BraidExecutor::set_spatial_coarsen_and_refine, "", "", "")
                        .add_method("set_refine", &T_BraidExecutor::set_refine, "", "", "")
                        .add_method("set_max_refinements", &T_BraidExecutor::set_max_refinements, "", "", "")
                        .add_method("set_access_level", &T_BraidExecutor::set_access_level, "", "", "")
                        .add_method("set_print_level", &T_BraidExecutor::set_print_level, "", "", "")
                        .add_method("set_print_file", &T_BraidExecutor::set_print_file, "", "", "")
                        .add_method("set_default_print_file", &T_BraidExecutor::set_default_print_file, "","", "")
                        .add_method("set_cycle_fmg", &T_BraidExecutor::set_cycle_fmg, "", "", "")
                        .add_method("set_cycle_nfmg", &T_BraidExecutor::set_cycle_nfmg, "", "", "")
                        .add_method("set_cycle_nfmgv", &T_BraidExecutor::set_cycle_nfmgv, "", "", "")
                        .add_method("set_cycle_type", &T_BraidExecutor::set_cycle_type, "", "", " ")
                        .add_method("set_sync", &T_BraidExecutor::set_sync, "", "", "")
                        .add_method("set_increase_max_levels", &T_BraidExecutor::set_increase_max_levels, "","", "")
                        .add_method("set_relax_only_cg", &T_BraidExecutor::set_relax_only_cg, "", "", "")
                        .add_method("set_agg_c_factor", &T_BraidExecutor::set_agg_c_factor, "", "", "")
                        .add_method("set_periodic", &T_BraidExecutor::set_periodic, "", "", "")
                        .add_method("set_final_fc_relax", &T_BraidExecutor::set_final_fc_relax, "", "", "")
                        .add_method("set_reverted_ranks", &T_BraidExecutor::set_reverted_ranks, "", "", "")
                        .add_method("set_richardson_estimation", &T_BraidExecutor::set_richardson_estimation, "", "", "")
                        .add_method("set_file_io_level", &T_BraidExecutor::set_file_io_level, "", "", "")
                        .add_method("set_c_relax_weight", &T_BraidExecutor::set_c_relax_weight, "", "", "")
                        .add_method("set_t_points_cutoff", &T_BraidExecutor::set_t_points_cutoff, "", "", "")
                        .add_method("set_full_residual_norm", &T_BraidExecutor::set_full_residual_norm, "","", "")
                        .add_method("set_time_grid", &T_BraidExecutor::set_time_grid, "", "", "")
                        .add_method("get_num_iteration", &T_BraidExecutor::get_num_iteration, "", "", "")
                        .add_method("get_c_factor", &T_BraidExecutor::get_c_factor, "", "", "")
                        .add_method("get_residual_norms", &T_BraidExecutor::get_residual_norms, "", "", "")
                        .add_method("get_num_level", &T_BraidExecutor::get_num_level, "", "", "")
                        .add_method("get_warm_restart", &T_BraidExecutor::get_warm_restart, "", "", "")
                        .add_method("get_distribution_lower", &T_BraidExecutor::get_distribution_lower, "","", "")
                        .add_method("get_distribution_upper", &T_BraidExecutor::get_distribution_upper, "","", "")
                        .add_method("get_id", &T_BraidExecutor::get_id, "", "", "")
                        .add_method("set_app", &T_BraidExecutor::set_app, "", "", "")
                        .add_method("set_paralog", &T_BraidExecutor::set_paralog, "", "", "")
                        .add_method("get_app", &T_BraidExecutor::get_app, "braid application", "", "")
                        .add_method("set_initializer", &T_BraidExecutor::set_initializer, "", "", "")
                        .add_method("set_norm_provider", &T_BraidExecutor::set_norm_provider, "", "", "")
                        .add_method("print_settings", &T_BraidExecutor::print_settings, "", "", "")
                        .add_method("print_summary", &T_BraidExecutor::print_summary, "", "", "")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name, "BraidExecutor", tag);
            }
        }


        template<typename TDomain>
        static void Domain(Registry &reg, string grp) {
            string suffix = GetDomainSuffix<TDomain>();
            string tag = GetDomainTag<TDomain>();
        }


        template<int dim>
        static void Dimension(Registry &reg, string grp) {
            string suffix = GetDimensionSuffix<dim>();
            string tag = GetDimensionTag<dim>();

        }


        template<typename TAlgebra>
        static void Algebra(Registry &reg, string grp) {
            string suffix = GetAlgebraSuffix<TAlgebra>();
            string tag = GetAlgebraTag<TAlgebra>();
        }


        static void Common(Registry &reg, string grp) {
            reg.add_function("get_virtual_memory_total", &get_virtual_memory_total, "", "", "");
            reg.add_function("get_virtual_memory_used", &get_virtual_memory_used, "", "", "");
            reg.add_function("get_virtual_memory_consumed", &get_virtual_memory_consumed, "", "", "");
            reg.add_function("get_physical_memory_total", &get_physical_memory_total, "", "", "");
            reg.add_function("get_physical_memory_used", &get_physical_memory_used, "", "", "");
            reg.add_function("get_physical_memory_consumed", &get_physical_memory_consumed, "", "","");
            reg.add_function("get_world_memory_consumed", &get_world_memory_consumed, "", "", "");
            reg.add_function("get_spatial_memory_consumed", &get_spatial_memory_consumed, "", "", "");
            reg.add_function("get_world_memory_distribution", &get_world_memory_distribution, "", "","");
            reg.add_function("get_spatial_memory_distribution", &get_spatial_memory_distribution, "","", "");
        }
    };
}// end namespace XBraidForUG4

extern "C" void
InitUGPlugin_XBraidForUG4(Registry *reg, string grp) {
    using namespace XBraidForUG4;
    grp.append("XBraidForUG4");
    // Space Time Communicator
    {
        typedef SpaceTimeCommunicator T_SpaceTimeCommunicator;
        string name = "SpaceTimeCommunicator";
        reg->add_class_<T_SpaceTimeCommunicator>(name, "XBraid", "")
                .add_constructor()
                .add_method("split", &T_SpaceTimeCommunicator::split)
                .add_method("unsplit", &T_SpaceTimeCommunicator::unsplit)
                .add_method("get_global_rank", &T_SpaceTimeCommunicator::get_global_rank)
                .add_method("get_spatial_rank", &T_SpaceTimeCommunicator::get_spatial_rank)
                .add_method("get_temporal_rank", &T_SpaceTimeCommunicator::get_temporal_rank)
                .add_method("get_global_size", &T_SpaceTimeCommunicator::get_global_size)
                .add_method("get_spatial_size", &T_SpaceTimeCommunicator::get_spatial_size)
                .add_method("get_temporal_size", &T_SpaceTimeCommunicator::get_temporal_size)
                .add_method("sleep", &T_SpaceTimeCommunicator::sleep)
                .set_construct_as_smart_pointer(true);
    }
    {// Paralog
        typedef ParallelLogger T_ParallelLogger;
        string name = "Paralog";
        reg->add_class_<T_ParallelLogger>(name, "XBraid", "")
                .add_constructor()
                .add_method("set_filename", &T_ParallelLogger::set_filename, "", "", "")
                .add_method("set_comm", &T_ParallelLogger::set_comm, "", "", "")
                .add_method("init", &T_ParallelLogger::init, "", "", "")
                .add_method("release", &T_ParallelLogger::release, "", "", "")
                .add_method("write", &T_ParallelLogger::write, "", "", "")
                .set_construct_as_smart_pointer(true);
    }
    { // ReplaceStandardStream
        typedef ReplaceStandardStream T_ReplaceStandardStream;
        string name = "ReplaceStandardStream";
        reg->add_class_<T_ReplaceStandardStream>(name, "XBraid", "")
                .add_constructor()
                .add_method("apply", &T_ReplaceStandardStream::apply, "", "", "")
                .add_method("undo", &T_ReplaceStandardStream::undo, "", "", "")
                .add_method("set_space_time_comm", &T_ReplaceStandardStream::set_space_time_comm, "", "", "")
                .set_construct_as_smart_pointer(true);
    }
    {
        typedef BraidTimer T_BraidTimer;
        reg->add_class_<T_BraidTimer>("BraidTimer", grp, "")
                .add_constructor()
                .add_method("start", &T_BraidTimer::start)
                .add_method("stop", &T_BraidTimer::stop)
                .add_method("get", &T_BraidTimer::get)
                .set_construct_as_smart_pointer(true);
    }
    try {
        RegisterCommon<Functionality>(*reg, grp);
        RegisterDimensionDependent<Functionality>(*reg, grp);
        RegisterDomainDependent<Functionality>(*reg, grp);
        RegisterAlgebraDependent<Functionality>(*reg, grp);
        RegisterDomainAlgebraDependent<Functionality>(*reg, grp);
    }
    UG_REGISTRY_CATCH_THROW(grp);
}}