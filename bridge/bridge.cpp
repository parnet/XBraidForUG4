#include "bridge/util.h"
#include "bridge/util_domain_algebra_dependent.h"

#include "core/braid_executor.hpp"
#include "core/space_time_communicator.hpp"
#include "core/replace_standard_stream.hpp"

#include "driver/gridfunction_base.hpp"
#include "driver/braid_theta_integrator_nonlinear.hpp"
#include "driver/braid_integrator.hpp"
#include "driver/braid_integrator_factory.hpp"
#include "driver/basic_driver.hpp"
#include "driver/braid_residual_stepper.hpp"

#include "factory/theta_integrator_factory.hpp"
#include "factory/fixed_step_theta_integrator_factory.hpp"
#include "factory/bdf_integrator_factory.hpp"

#include "initializer/start_value_initializer.hpp"
#include "initializer/zero_Initializer.hpp"
#include "initializer/random_value_initializer.hpp"

#include "integrator/bdf_integrator.hpp"
#include "integrator/theta_conststep_integrator.hpp"
#include "integrator/theta_conststep_integrator_nl.hpp"
#include "integrator/theta_integrator_nl.hpp"
#include "integrator/bdf_integrator_nl.hpp"
#include "integrator/theta_single_timestep.hpp"
#include "integrator/experimental_timestep.hpp"

#include "interface/initializer.hpp"
#include "interface/residual_timeintegrator.hpp"
#include "interface/integrator_factory.hpp"
#include "interface/spatial_norm.hpp"

#include "factory/limex_integrator_factory.hpp"
#include "factory/simple_integrator_factory.hpp"
#include "factory/const_step_linear_time_integrator_factory.hpp"
#include "factory/linear_time_integrator_factory.hpp"

#include "spatial_norm/euclidian_norm.hpp"

#include "util/braid_timer.hpp"
#include "util/io_gridfunction.hpp"
#include "util/memory_observer.hpp"
#include "util/parallel_io_gridfunction.hpp"
#include "util/parallel_logger.hpp"
#include "util/world_memory.hpp"

#include "observer/collector_observer.hpp"
#include "observer/io_observer.hpp"
#include "observer/io_process_observer.hpp"
#include "observer/matlab_observer.hpp"
//#include "observer/mock_observer.hpp"
#include "core/discard_standard_stream.hpp"
#include "observer/userdata_process_evaluate_observer.hpp"
#include "observer/vtk_observer.hpp"
#include "observer/vtk_process_observer.hpp"
#include "observer/xb_collector_observer.hpp"

#include "bridge/bridge.h"


#ifdef XBraidPoroelasticity
#include "problemset/poroelasticity/bridge.hpp"
#endif

namespace ug {

    namespace xbraid {

    struct Functionality {
        template<typename TDomain, typename TAlgebra>
        static void DomainAlgebra(Registry &reg, std::string grp) {


            std::string suffix = bridge::GetDomainAlgebraSuffix<TDomain, TAlgebra>();
            std::string tag = bridge::GetDomainAlgebraTag<TDomain, TAlgebra>();
            // SpatialGridTransfer
            {
                using T_SpatialGridTransfer = SpatialGridTransfer<TDomain, TAlgebra>;
                std::string name = std::string("SpatialGridTransfer").append(suffix);
                reg.add_class_<T_SpatialGridTransfer>(name,grp)
                        .add_constructor()
                        .add_method("set_transfer", &T_SpatialGridTransfer::set_transfer, "", "", "")
                        .add_method("make_nontop", &T_SpatialGridTransfer::make_nontop, "", "", "")
                        .add_method("prolongate", &T_SpatialGridTransfer::prolongate, "", "", "")
                        .add_method("restrict", &T_SpatialGridTransfer::restrict, "", "", "")
                        .add_method("set_domain", &T_SpatialGridTransfer::set_domain, "", "", "")
                        .add_method("set_approx_space", &T_SpatialGridTransfer::set_approx_space, "", "", "")
                        .add_method("set_prolongation", &T_SpatialGridTransfer::set_prolongation, "", "", "")
                        .add_method("set_restriction", &T_SpatialGridTransfer::set_restriction, "", "", "")
                        .add_method("init", &T_SpatialGridTransfer::init, "", "", "")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name, "SpatialGridTransfer", tag);
            }




            /***********************************************************************************************************
             *  Integrator Factories
             **********************************************************************************************************/
            {
                // Interface
                {
                    using T_IntegratorFactory = IntegratorFactory<TDomain, TAlgebra> ;
                    std::string name = std::string("IntegratorFactory").append(suffix);
                    reg.add_class_<T_IntegratorFactory>(name, grp);
                    reg.add_class_to_group(name, "IntegratorFactory", tag);
                }
                // Integrator Factory
                {
                    using T_LimexFactory = LimexFactory<TDomain, TAlgebra> ;
                    using T_IntegratorFactory = IntegratorFactory<TDomain, TAlgebra> ;

                    std::string name = std::string("LimexFactory").append(suffix);
                    reg.add_class_<T_LimexFactory, T_IntegratorFactory>(name,grp)
                            .add_constructor()
                            .add_method("set_domain", &T_LimexFactory::set_domain, "", "", "")
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
                // LinearTimeIntegratorFactory
                {
                    using T_LinearTimeIntegratorFactory = LinearTimeIntegratorFactory<TDomain, TAlgebra> ;
                    using T_IntegratorFactory = IntegratorFactory<TDomain, TAlgebra> ;
                    std::string name = std::string("LinearTimeIntegratorFactory").append(suffix);
                    reg.add_class_<T_LinearTimeIntegratorFactory, T_IntegratorFactory>(name,grp)
                            .add_constructor()
                            .add_method("set_time_disc", &T_LinearTimeIntegratorFactory::set_time_disc, "", "", "")
                            .add_method("set_solver", &T_LinearTimeIntegratorFactory::set_solver, "", "", "")
                            .add_method("create_time_integrator", &T_LinearTimeIntegratorFactory::create_time_integrator,"", "", "")
                            .add_method("create_level_time_integrator",&T_LinearTimeIntegratorFactory::create_level_time_integrator, "", "", "")
                            .set_construct_as_smart_pointer(true);
                    reg.add_class_to_group(name, "LinearTimeIntegratorFactory", tag);
                }
                // ConstStepLinearTimeIntegratorFactory
                {
                    using T_ConstStepLinearTimeIntegratorFactory = ConstStepLinearTimeIntegratorFactory<TDomain, TAlgebra> ;
                    using T_IntegratorFactory = IntegratorFactory<TDomain, TAlgebra> ;
                    std::string name = std::string("ConstStepLinearTimeIntegratorFactory").append(suffix);
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
                // ThetaIntegratorFactory
                {
                    using ThetaIntegratorFactory = ThetaIntegratorFactory<TDomain, TAlgebra> ;
                    using T_IntegratorFactory = IntegratorFactory<TDomain, TAlgebra> ;
                    std::string name = std::string("ThetaIntegratorFactory").append(suffix);
                    reg.add_class_<ThetaIntegratorFactory, T_IntegratorFactory>(name, grp)
                            .add_constructor()
                            .add_method("set_domain", &ThetaIntegratorFactory::set_domain, "", "","")
                            .add_method("set_solver", &ThetaIntegratorFactory::set_solver, "", "", "")
                            .add_method("set_theta", &ThetaIntegratorFactory::set_theta, "", "", "")
                            .add_method("set_level_theta", &ThetaIntegratorFactory::set_level_theta, "", "", "")
                            .add_method("create_time_integrator", &ThetaIntegratorFactory::create_time_integrator, "", "","")
                            .add_method("create_level_time_integrator",&ThetaIntegratorFactory::create_level_time_integrator, "", "", "")
                            .set_construct_as_smart_pointer(true);
                    reg.add_class_to_group(name, "ThetaIntegratorFactory", tag);
                }
                // FixedStepThetaIntegratorFactory
                {
                    using T_FixedStepThetaIntegratorFactory = FixedStepThetaIntegratorFactory<TDomain, TAlgebra> ;
                    using TBase = IntegratorFactory<TDomain, TAlgebra> ;
                    std::string name = std::string("FixedStepThetaIntegratorFactory").append(suffix);
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
                // BDF_IntegratorFactory
                {
                    using T_BDF_IntegratorFactory = BDF_IntegratorFactory<TDomain, TAlgebra> ;
                    using TBase = IntegratorFactory<TDomain, TAlgebra> ;
                    std::string name = std::string("BDF_IntegratorFactory").append(suffix);
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
                // SimpleIntegratorFactory
                {
                    using T_SimpleIntegratorFactory = SimpleIntegratorFactory<TDomain, TAlgebra> ;
                    using T_IntegratorFactory = IntegratorFactory<TDomain, TAlgebra> ;
                    std::string name = std::string("SimpleIntegratorFactory").append(suffix);
                    reg.add_class_<T_SimpleIntegratorFactory, T_IntegratorFactory>(name,grp)
                            .add_constructor()
                            .add_method("set_domain", &T_SimpleIntegratorFactory::set_domain, "", "", "")
                            .add_method("set_solver", &T_SimpleIntegratorFactory::set_solver, "", "", "")
                            .add_method("set_dt_min", &T_SimpleIntegratorFactory::set_dt_min, "", "end time", "")
                            .add_method("set_dt_max", &T_SimpleIntegratorFactory::set_dt_max, "", "", "")
                            .add_method("create_time_integrator", &T_SimpleIntegratorFactory::create_time_integrator,"", "", "")
                            .add_method("create_level_time_integrator",&T_SimpleIntegratorFactory::create_level_time_integrator, "", "", "")
                            .set_construct_as_smart_pointer(true);
                    reg.add_class_to_group(name, "SimpleIntegratorFactory", tag);
                }
            }



            /***********************************************************************************************************
             *  Integrator Methods
             **********************************************************************************************************/
            {
                // Interface
                {
                    using  T_IResidualTimeIntegrator = IResidualTimeIntegrator<TDomain, TAlgebra> ;
                    std::string name = std::string("IResidualTimeIntegrator").append(suffix);
                    reg.add_class_<T_IResidualTimeIntegrator>(name, grp);
                    reg.add_class_to_group(name, "IResidualTimeIntegrator", tag);
                }
                // BDF_Integrator
                {
                    using T_BDF_Integrator = BDF_Integrator<TDomain, TAlgebra> ;
                    using TBase = ITimeIntegrator<TDomain, TAlgebra> ;
                    std::string name = std::string("BDF_Integrator").append(suffix);
                    reg.add_class_<T_BDF_Integrator, TBase>(name, grp)
                            .add_constructor()
                            .add_method("set_domain", &T_BDF_Integrator::set_domain, "", "", "")
                            .add_method("set_solver", &T_BDF_Integrator::set_solver, "", "", "")
                            .add_method("set_order", &T_BDF_Integrator::set_order, "", "", "")
                            .add_method("set_reassemble_threshold", &T_BDF_Integrator::set_reassemble_threshold, "","", "")
                            .set_construct_as_smart_pointer(true);
                    reg.add_class_to_group(name, "BDF_Integrator", tag);
                }
                // BDF_IntegratorNL
                {
                    using T_NLBDFIntegrator = BDF_IntegratorNL<TDomain, TAlgebra> ;
                    using T_ITimeIntegrator = ITimeIntegrator<TDomain, TAlgebra> ;
                    std::string name = std::string("BDF_IntegratorNL").append(suffix);
                    reg.add_class_<T_NLBDFIntegrator, T_ITimeIntegrator>(name, grp)
                            .add_constructor()
                            .add_method("set_domain", &T_NLBDFIntegrator::set_domain, "", "", "")
                            .add_method("set_solver", &T_NLBDFIntegrator::set_solver, "", "", "")
                            .add_method("set_order", &T_NLBDFIntegrator::set_order, "", "", "")
                            .set_construct_as_smart_pointer(true);
                    reg.add_class_to_group(name, "BDF_IntegratorNL", tag);
                }
                // ThetaConstStepIntegrator
                {
                using T_ThetaConstStepIntegrator = ThetaConstStepIntegrator<TDomain, TAlgebra> ;
                using T_Base = ITimeIntegrator<TDomain, TAlgebra> ;
                std::string name = std::string("ThetaConstStepIntegrator").append(suffix);
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
                // ThetaIntegratorNL
                {
                    using T_ThetaIntegratorNL = ThetaIntegratorNL<TDomain, TAlgebra> ;
                    using T_ITimeIntegrator = ITimeIntegrator<TDomain, TAlgebra> ;
                    std::string name = std::string("ThetaIntegratorNL").append(suffix);
                    reg.add_class_<T_ThetaIntegratorNL, T_ITimeIntegrator>(name, grp)
                            .add_constructor()
                            .add_method("set_domain", &T_ThetaIntegratorNL::set_domain, "", "", "")
                            .add_method("set_solver", &T_ThetaIntegratorNL::set_solver, "", "", "")
                            .add_method("set_theta", &T_ThetaIntegratorNL::set_theta, "", "", "")
                            .add_method("apply", &T_ThetaIntegratorNL::apply, "", "", "")
                            .set_construct_as_smart_pointer(true);
                    reg.add_class_to_group(name, "ThetaIntegratorNL", tag);
                }
                // ThetaConstStepIntegratorNL
                {
                    using T_ThetaConstStepIntegratorNL = ThetaConstStepIntegratorNL<TDomain, TAlgebra> ;
                    using T_ITimeIntegrator = ITimeIntegrator<TDomain, TAlgebra> ;
                    std::string name = std::string("ThetaConstStepIntegratorNL").append(suffix);
                    reg.add_class_<T_ThetaConstStepIntegratorNL, T_ITimeIntegrator>(name, grp)
                            .add_constructor()
                            .add_method("set_domain", &T_ThetaConstStepIntegratorNL::set_domain, "", "", "")
                            .add_method("set_solver", &T_ThetaConstStepIntegratorNL::set_solver, "", "", "")
                            .add_method("set_theta", &T_ThetaConstStepIntegratorNL::set_theta, "", "", "")
                            .add_method("set_num_steps", &T_ThetaConstStepIntegratorNL::set_num_steps, "", "", "")
                            .set_construct_as_smart_pointer(true);
                    reg.add_class_to_group(name, "ThetaConstStepIntegratorNL", tag);
                }
                // ThetaSingleTimeStep
                {
                    using  T_GridFunction = GridFunction<TDomain, TAlgebra> ;
                    using SP_GridFunction =  SmartPtr<T_GridFunction> ;
                    using CP_GridFunction = ConstSmartPtr<T_GridFunction> ;
                    using T_ThetaSingleTimeStep = ThetaSingleTimeStep<TDomain, TAlgebra> ;
                    using T_IResidualTimeIntegrator = IResidualTimeIntegrator<TDomain, TAlgebra> ;
                    std::string name = std::string("ThetaSingleTimeStep").append(suffix);
                    reg.add_class_<T_ThetaSingleTimeStep, T_IResidualTimeIntegrator>(name, grp)
                            .add_constructor()
                            .add_method("set_domain", &T_ThetaSingleTimeStep::set_domain, "", "", "")
                            .add_method("set_solver", &T_ThetaSingleTimeStep::set_solver, "", "", "")
                            .add_method("set_theta", &T_ThetaSingleTimeStep::set_theta, "", "", "")
                            .add_method("apply", static_cast<bool (T_ThetaSingleTimeStep::*)(SP_GridFunction , double , CP_GridFunction , double)>(&T_ThetaSingleTimeStep::apply) , "", "", "")
                            //.add_method("apply", bool (&T_ThetaSingleTimeStep::apply)(SP_GridFunction , double , CP_GridFunction , double )) , "", "", "")
                            .add_method("set_reassemble_threshold", &T_ThetaSingleTimeStep::set_reassemble_threshold, "","", "")
                            .set_construct_as_smart_pointer(true);
                    reg.add_class_to_group(name, "ThetaSingleTimeStep", tag);
                }
                // ø Braid Time Stepper
                {
                   /* using T_BraidTimeStepper = ThetaSingleTimeStep<TDomain, TAlgebra> ;
                    using T_BraidGridFunctionBase = IResidualTimeIntegrator<TDomain, TAlgebra> ;
                    std::string name = std::string("BraidResidualStepper").append(suffix);
                    reg.add_class_<T_BraidTimeStepper, T_BraidGridFunctionBase>(name, grp)
                            .add_constructor()
                            //.add_method("print_settings", &T_BraidTimeStepper::print_settings, "", "", "")
                            //.add_method("set_approx_space", &T_BraidTimeStepper::set_approx_space, "", "", "")
                            //.add_method("set_adapt_convergence", &T_BraidTimeStepper::setAdaptConv, "", "")
                            .add_method("set_domain", &T_BraidTimeStepper::set_domain, "", "")
                            .add_method("set_solver", &T_BraidTimeStepper::set_solver, "", "")
                            //.add_method("set_force_convergence", &T_BraidTimeStepper::setForceConv, "", "", "")
                            .set_construct_as_smart_pointer(true);
                    reg.add_class_to_group(name, "BraidResidualStepper", tag);*/
                }
                // ExperimentalTimeStep
                {
                    using  T_GridFunction = GridFunction<TDomain, TAlgebra> ;
                    using SP_GridFunction =  SmartPtr<T_GridFunction> ;
                    using CP_GridFunction = ConstSmartPtr<T_GridFunction> ;
                    using T_ExperimentalTimeStep = ExperimentalTimeStep<TDomain, TAlgebra> ;
                    using T_IResidualTimeIntegrator = IResidualTimeIntegrator<TDomain, TAlgebra> ;
                    std::string name = std::string("ExperimentalTimeStep").append(suffix);
                    reg.add_class_<T_ExperimentalTimeStep, T_IResidualTimeIntegrator>(name, grp)
                            .add_constructor()
                            .add_method("set_domain", &T_ExperimentalTimeStep::set_domain, "", "", "")
                            .add_method("set_solver", &T_ExperimentalTimeStep::set_solver, "", "", "")
                            .add_method("set_theta", &T_ExperimentalTimeStep::set_theta, "", "", "")
                            .add_method("set_rhs", &T_ExperimentalTimeStep::set_rhs, "", "", "")
                            .add_method("set_jacobian", &T_ExperimentalTimeStep::set_jacobian, "", "", "")
                            .add_method("apply", static_cast<bool (T_ExperimentalTimeStep::*)(SP_GridFunction , double , CP_GridFunction , double)>(&T_ExperimentalTimeStep::apply) , "", "", "")
                            //.add_method("apply", bool (&T_ThetaSingleTimeStep::apply)(SP_GridFunction , double , CP_GridFunction , double )) , "", "", "")
                            .add_method("set_reassemble_threshold", &T_ExperimentalTimeStep::set_reassemble_threshold, "","", "")
                            .set_construct_as_smart_pointer(true);
                    reg.add_class_to_group(name, "ExperimentalTimeStep", tag);
                }
            }



            /***********************************************************************************************************
             *  Utility
             **********************************************************************************************************/
            {
                // IOGridFunction
                {
                    using T_IOGridFunction = IOGridFunction<TDomain, TAlgebra> ;
                    std::string name = std::string("IOGridFunction").append(suffix);
                    reg.add_class_<T_IOGridFunction>(name, grp)
                            .add_constructor()
                            .add_method("write", &T_IOGridFunction::write, "", "", "")
                            .add_method("read", &T_IOGridFunction::read, "", "", "")
                            .set_construct_as_smart_pointer(true);
                    reg.add_class_to_group(name, "IOGridFunction", tag);
                }
                // PIOGridFunction
                {
                    using T_PIOGridFunction = PIOGridFunction<TDomain, TAlgebra> ;
                    std::string name = std::string("PIOGridFunction").append(suffix);
                    reg.add_class_<T_PIOGridFunction>(name, grp)
                            .add_constructor()
                            .add_method("write", &T_PIOGridFunction::write, "", "", "")
                            .add_method("read", &T_PIOGridFunction::read, "", "", "")
                            .set_construct_as_smart_pointer(true);
                    reg.add_class_to_group(name, "PIOGridFunction", tag);
                }
            }



            /***********************************************************************************************************
             *  Initializer
             **********************************************************************************************************/
            {
                // Interface
                {
                    using T_BraidInitializer = BraidInitializer<TDomain, TAlgebra> ;
                    std::string name = std::string("BraidInitializer").append(suffix);
                    reg.add_class_<T_BraidInitializer>(name, grp)
                        .add_method("set_start_values", &T_BraidInitializer::set_start_values, "", "", "")
                        .set_construct_as_smart_pointer(true);
                    reg.add_class_to_group(name, "BraidInitializer", tag);
                }
                // StartValueInitializer
                {
                    using T_StartValueInitializer = GridFunctionInitializer<TDomain, TAlgebra> ;
                    using T_BraidInitializer = BraidInitializer<TDomain, TAlgebra> ;
                    std::string name = std::string("GridFunctionInitializer").append(suffix);
                    reg.add_class_<T_StartValueInitializer, T_BraidInitializer>(name, grp)
                            .add_constructor()
                            .set_construct_as_smart_pointer(true);
                    reg.add_class_to_group(name, "GridFunctionInitializer", tag);
                }
                // ZeroInitializer
                {
                    using T_ZeroValueInitializer = ZeroInitializer<TDomain, TAlgebra> ;
                    using T_BraidInitializer = BraidInitializer<TDomain, TAlgebra> ;
                    std::string name = std::string("ZeroInitializer").append(suffix);
                    reg.add_class_<T_ZeroValueInitializer, T_BraidInitializer>(name, grp)
                            .add_constructor()
                            .set_construct_as_smart_pointer(true);
                    reg.add_class_to_group(name, "ZeroInitializer", tag);
                }
                // RandomValueInitializer
                {
                    using T_RandomValueInitializer = RandomValueInitializer<TDomain, TAlgebra> ;
                    using T_BraidInitializer = BraidInitializer<TDomain, TAlgebra> ;
                    std::string name = std::string("RandomValueInitializer").append(suffix);
                    reg.add_class_<T_RandomValueInitializer, T_BraidInitializer>(name, grp)
                            .add_constructor()
                            .add_method("set_parameter_uniform", &T_RandomValueInitializer::set_parameter_uniform, "", "","")
                            .add_method("set_parameter_normal", &T_RandomValueInitializer::set_parameter_normal, "", "","")
                            //2025-04 .add_method("set_domain", &T_RandomValueInitializer::set_domain, "", "","")
                            .set_construct_as_smart_pointer(true);
                    reg.add_class_to_group(name, "RandomValueInitializer", tag);
                }
            }

            /***********************************************************************************************************
             *  Observer
             **********************************************************************************************************/
            {
                // Interface
                {
                    using T_IXBraidTimeIntegratorObserver = IXBraidTimeIntegratorObserver<TDomain, TAlgebra> ;
                    using T_ITimeIntegratorObserver = ITimeIntegratorObserver<TDomain, TAlgebra> ;
                    std::string name = std::string("XBraidTimeIntegratorObserver").append(suffix);
                    reg.add_class_<T_IXBraidTimeIntegratorObserver, T_ITimeIntegratorObserver>(name, grp)
                            .add_method("write", &T_IXBraidTimeIntegratorObserver::write, "", "","")
                            .set_construct_as_smart_pointer(true);
                    reg.add_class_to_group(name, "XBraidTimeIntegratorObserver", tag);
                }

                // TimeIntegratorObserverCollector (may be deprecated)
                {
                    using T_TimeIntegratorObserverCollector = TimeIntegratorObserverCollector<TDomain, TAlgebra> ;
                    using T_ITimeIntegratorObserver = ITimeIntegratorObserver<TDomain, TAlgebra> ;
                    std::string name_multi = std::string("TimeIntegratorObserverCollector").append(suffix);
                    reg.add_class_<T_TimeIntegratorObserverCollector, T_ITimeIntegratorObserver>(name_multi, grp)
                            .add_constructor()
                            .add_method("attach_observer", &T_TimeIntegratorObserverCollector::attach_observer, "", "", "")
                            .set_construct_as_smart_pointer(true);
                    reg.add_class_to_group(name_multi, "TimeIntegratorObserverCollector", tag);
                }

                // XBraid Observer Collector
                {
                    using T_XBraidTimeIntegratorObserverCollector = XBraidTimeIntegratorObserverCollector<TDomain, TAlgebra> ;
                    using T_IXBraidTimeIntegratorObserver = IXBraidTimeIntegratorObserver<TDomain, TAlgebra> ;
                    using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
                    using SP_GridFunction = SmartPtr<T_GridFunction> ;
                    std::string name_multi = std::string("XBraidTimeIntegratorObserverCollector").append(suffix);
                    reg.add_class_<T_XBraidTimeIntegratorObserverCollector, T_IXBraidTimeIntegratorObserver>(name_multi, grp)
                            .add_constructor()
                            .add_method("attach_observer", &T_XBraidTimeIntegratorObserverCollector::attach_observer, "", "", "")
                            .add_method("attach_common_observer", &T_XBraidTimeIntegratorObserverCollector::attach_common_observer, "", "", "")
                            .add_method("step_process",(bool (T_XBraidTimeIntegratorObserverCollector::*)(SP_GridFunction u, int, number, number)) &T_XBraidTimeIntegratorObserverCollector::step_process, "", "", "")
                            .add_method("step_process", (bool (T_XBraidTimeIntegratorObserverCollector::*)(SP_GridFunction u, int, number, number, int, int)) &T_XBraidTimeIntegratorObserverCollector::step_process, "","", "")
                            .set_construct_as_smart_pointer(true);
                    reg.add_class_to_group(name_multi, "XBraidTimeIntegratorObserverCollector", tag);
                }

                // VTK_Observer
                {
                    using T_VTK_Observer = VTK_Observer<TDomain, TAlgebra> ;
                    using T_ITimeIntegratorObserver = ITimeIntegratorObserver<TDomain, TAlgebra> ;
                    using T_VTKOutput = VTKOutput<TDomain::dim> ;
                    using SP_VTKOutput = SmartPtr<T_VTKOutput> ;
                    std::string name_multi = std::string("VTK_Observer").append(suffix);
                    reg.add_class_<T_VTK_Observer, T_ITimeIntegratorObserver>(name_multi, grp)
                            .template add_constructor<void (*)(SP_VTKOutput, const char *)>("")
                            .add_method("step_process", &T_VTK_Observer::step_process, "", "", "")
                            .add_method("write_time_pvd", &T_VTK_Observer::write_time_pvd, "", "", "")
                            .add_method("set_filename", &T_VTK_Observer::set_filename, "", "", "")
                            .set_construct_as_smart_pointer(true);
                    reg.add_class_to_group(name_multi, "VTK_Observer", tag);
                }

                // MATLAB_Observer
                {
                    using T_MATLAB_Observer = MATLAB_Observer<TDomain, TAlgebra> ;
                    using T_ITimeIntegratorObserver = ITimeIntegratorObserver<TDomain, TAlgebra> ;

                    using SP_ParallelLogger = SmartPtr<ParallelLogger> ;
                    std::string name_multi = std::string("MATLAB_Observer").append(suffix);
                    reg.add_class_<T_MATLAB_Observer, T_ITimeIntegratorObserver>(name_multi, grp)
                    .template add_constructor<void (*)(SP_ParallelLogger)>(", ")
                            .set_construct_as_smart_pointer(true);
                    reg.add_class_to_group(name_multi, "MATLAB_Observer", tag);
                }

                // XBraid VTK Process Observer
                {
                    using T_VTK_ProcessObserver = VTK_ProcessObserver<TDomain, TAlgebra> ;
                    using T_IXBraidTimeIntegratorObserver = IXBraidTimeIntegratorObserver<TDomain, TAlgebra> ;

                    using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
                    using SP_GridFunction = SmartPtr<T_GridFunction> ;

                    using T_VTKOutput = VTKOutput<TDomain::dim> ;
                    using SP_VTKOutput = SmartPtr<T_VTKOutput> ;

                    std::string name_multi = std::string("VTK_ProcessObserver").append(suffix);
                    reg.add_class_<T_VTK_ProcessObserver, T_IXBraidTimeIntegratorObserver>(name_multi, grp)
                            .template add_constructor<void (*)(SP_VTKOutput, const char *)>("")
                            .add_method("step_process", (bool (T_VTK_ProcessObserver::*)(SP_GridFunction u, int,number,number) ) &T_VTK_ProcessObserver::step_process, "","", "")
                            .add_method("step_process", (bool (T_VTK_ProcessObserver::*)(SP_GridFunction u, int,number,number, int,int) ) &T_VTK_ProcessObserver::step_process, "","", "")
                            .set_construct_as_smart_pointer(true);
                    reg.add_class_to_group(name_multi, "VTK_ProcessObserver", tag);
                }

                // EvalObserver / UserdataProcessEvaluateObserver
                {
                    using T_UserdataProcessEvaluateObserver = UserdataProcessEvaluateObserver<TDomain, TAlgebra> ;
                    using T_IXBraidTimeIntegratorObserver = IXBraidTimeIntegratorObserver<TDomain, TAlgebra> ;
                    std::string name_eval = std::string("EvalObserver").append(suffix);
                    reg.add_class_<T_UserdataProcessEvaluateObserver, T_IXBraidTimeIntegratorObserver>(name_eval, grp)
                            .add_constructor()
                            .add_method("set_filename",static_cast<void (T_UserdataProcessEvaluateObserver::*)(const char *)>(&T_UserdataProcessEvaluateObserver::set_filename),"","", "")
                            .add_method("set_generator_component", &T_UserdataProcessEvaluateObserver::set_generator_component, "","", "")
                            .add_method("set_vector_generator", static_cast<void (T_UserdataProcessEvaluateObserver::*)(const char *)>(&T_UserdataProcessEvaluateObserver::set_vector_generator), "", "","")
                            .add_method("set_domain", &T_UserdataProcessEvaluateObserver::set_domain, "", "","")
                            .add_method("set_relative", &T_UserdataProcessEvaluateObserver::set_relative, "","", "")
                            .add_method("write_time_pvd", &T_UserdataProcessEvaluateObserver::write_time_pvd, "", "","")
                            .set_construct_as_smart_pointer(true);
                    reg.add_class_to_group(name_eval, "EvalObserver", tag);
                }
            }
            /***********************************************************************************************************
             *  Spatial Norm
             **********************************************************************************************************/
            {
                //Interface
                {
                    using T_BraidSpatialNorm = BraidSpatialNorm<TDomain, TAlgebra> ;
                    std::string name = std::string("BraidSpatialNorm").append(suffix);
                    reg.add_class_<T_BraidSpatialNorm>(name, grp);
                    reg.add_class_to_group(name, "BraidSpatialNorm", tag);
                }
                // Euclidian Norm / l2
                {
                    using T_BraidEuclidianNorm = BraidEuclidianNorm<TDomain, TAlgebra> ;
                    using T_BraidSpatialNorm = BraidSpatialNorm<TDomain, TAlgebra> ;
                    std::string name = std::string("BraidEuclidianNorm").append(suffix);
                    reg.add_class_<T_BraidEuclidianNorm, T_BraidSpatialNorm>(name, grp)
                            .add_constructor()
                            .add_method("norm", &T_BraidEuclidianNorm::norm, "", "", "")
                            .set_construct_as_smart_pointer(true);;
                    reg.add_class_to_group(name, "BraidEuclidianNorm", tag);
                }
            }



            // Braid Grid Function Base
            {
                using T_BraidGridFunctionBase = BraidGridFunctionBase<TDomain, TAlgebra> ;
                std::string name = std::string("BraidGridFunctionBase").append(suffix);
                reg.add_class_<T_BraidGridFunctionBase>(name, grp)
                        .add_method("init", &T_BraidGridFunctionBase::init, "", "", "")
                        .add_method("set_start_time", &T_BraidGridFunctionBase::set_start_time, "", "")
                        .add_method("set_end_time", &T_BraidGridFunctionBase::set_end_time, "", "", "")
                        .add_method("set_number_of_timesteps", &T_BraidGridFunctionBase::set_number_of_timesteps, "", "", "")
                        .add_method("set_time_values", static_cast<void (T_BraidGridFunctionBase::*)(double, double, int)>(&T_BraidGridFunctionBase::set_time_values), "", "", "")
                        .add_method("set_start_vector", &T_BraidGridFunctionBase::set_start_vector, "", "", "")
                        .add_method("set_norm_provider", &T_BraidGridFunctionBase::set_norm_provider, "", "", "")
                        .add_method("attach_xbraid_observer", &T_BraidGridFunctionBase::attach_xbraid_observer, "", "", "") // todo rename to process observer
                        .add_method("attach_observer", &T_BraidGridFunctionBase::attach_observer, "", "", "")
                        .add_method("set_max_levels", &T_BraidGridFunctionBase::set_max_levels, "", "")
                        .add_method("set_domain", &T_BraidGridFunctionBase::set_domain, "", "")
                        .add_method("set_initializer", &T_BraidGridFunctionBase::set_initializer, "", "");
                reg.add_class_to_group(name, "BraidGridFunctionBase", tag);
            }
            // Braid Time Integrator
            {
                using T_BraidIntegrator = BraidIntegrator<TDomain, TAlgebra> ;
                using T_BraidGridFunctionBase = BraidGridFunctionBase<TDomain, TAlgebra> ;
                std::string name = std::string("BraidIntegrator").append(suffix);
                reg.add_class_<T_BraidIntegrator, T_BraidGridFunctionBase>(name, grp)
                        .add_constructor()
                        .add_method("print_settings", &T_BraidIntegrator::print_settings, "", "", "")
                        .add_method("set_ref_factor", &T_BraidIntegrator::set_ref_factor, "", "", "")
                        .add_method("set_threshold", &T_BraidIntegrator::set_threshold, "", "", "")
                        .add_method("set_integrator", &T_BraidIntegrator::set_integrator, "", "", "")
                        .add_method("set_default_integrator", &T_BraidIntegrator::set_default_integrator, "", "", "")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name, "BraidIntegrator", tag);
            }
            // BraidNLIntegrator
            {
                using T_BraidNLIntegrator = BraidNLIntegrator<TDomain, TAlgebra> ;
                using T_BraidGridFunctionBase = BraidGridFunctionBase<TDomain, TAlgebra> ;
                std::string name = std::string("BraidNLIntegrator").append(suffix);
                reg.add_class_<T_BraidNLIntegrator, T_BraidGridFunctionBase>(name, grp)
                        .add_constructor()
                        .add_method("print_settings", &T_BraidNLIntegrator::print_settings, "", "", "")

                        .add_method("set_threshold", &T_BraidNLIntegrator::set_threshold, "", "", "")
                        .add_method("set_integrator", &T_BraidNLIntegrator::set_integrator, "", "", "")
                        .add_method("set_default_integrator", &T_BraidNLIntegrator::set_default_integrator, "", "", "")
                        .add_method("set_conv_check", &T_BraidNLIntegrator::set_conv_check, "", "", "")
                        .add_method("set_tol", &T_BraidNLIntegrator::set_tol, "", "", "")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name, "BraidNLIntegrator", tag);
            }
            // BraidIntegratorFactory
            {
                using T_BraidIntegratorFactory = BraidIntegratorFactory<TDomain, TAlgebra> ;
                using T_BraidGridFunctionBase = BraidGridFunctionBase<TDomain, TAlgebra> ;
                std::string name = std::string("BraidIntegratorFactory").append(suffix);
                reg.add_class_<T_BraidIntegratorFactory, T_BraidGridFunctionBase>(name, grp)
                        .add_constructor()
                        .add_method("print_settings", &T_BraidIntegratorFactory::print_settings, "", "", "")
                        .add_method("set_default_integrator", &T_BraidIntegratorFactory::set_default_integrator, "", "", "")
                        .add_method("set_integrator", &T_BraidIntegratorFactory::set_integrator, "", "", "")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name, "BraidIntegratorFactory", tag);
            }

            // Basic Driver
            {
                using T_BasicDriver = BasicDriver<TDomain, TAlgebra> ;
                using T_BraidGridFunctionBase = BraidGridFunctionBase<TDomain, TAlgebra> ;
                std::string name = std::string("BasicDriver").append(suffix);
                reg.add_class_<T_BasicDriver, T_BraidGridFunctionBase>(name, grp)
                        .add_constructor()
                        .add_method("print_settings", &T_BasicDriver::print_settings, "", "", "")
                        .add_method("set_domain", &T_BasicDriver::set_domain, "", "", "")
                        .add_method("set_integrator", &T_BasicDriver::set_integrator, "", "", "")
                        .add_method("set_default_integrator", &T_BasicDriver::set_default_integrator, "", "","")
#ifdef FEATURE_SPATIAL_REFINE
                        .add_method("set_spatial_grid_transfer", &T_BasicDriver::set_spatial_grid_transfer, "", "", "")
                        .add_method("set_level_num_ref", &T_BasicDriver::set_level_num_ref, "", "", "")
#endif


                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name, "BasicDriver", tag);
            }
            // BraidExecutor
            {
                using T_BraidExecutor = BraidExecutor<TDomain, TAlgebra> ;
                using SP_SpaceTimeCommunicator = SmartPtr<SpaceTimeCommunicator> ;
                using T_BraidGridFunctionBase = BraidGridFunctionBase<TDomain, TAlgebra> ;
                using SP_BraidGridFunctionBase = SmartPtr<T_BraidGridFunctionBase> ;
                using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
                using SP_GridFunction = SmartPtr<T_GridFunction> ;
                std::string name = std::string("BraidExecutor").append(suffix);
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
                        .add_method("set_driver", &T_BraidExecutor::set_driver, "", "", "")
                        .add_method("set_parallel_logger", &T_BraidExecutor::set_parallel_logger, "", "", "")
                        .add_method("get_driver", &T_BraidExecutor::get_driver, "braid driver", "", "")
                        .add_method("set_initializer", &T_BraidExecutor::set_initializer, "", "", "")
                        .add_method("set_norm_provider", &T_BraidExecutor::set_norm_provider, "", "", "")
                        .add_method("print_settings", &T_BraidExecutor::print_settings, "", "", "")
                        .add_method("print_summary", &T_BraidExecutor::print_summary, "", "", "")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name, "BraidExecutor", tag);
            }
        }


        template<typename TDomain>
        static void Domain(Registry &reg, std::string grp) {
            std::string suffix = bridge::GetDomainSuffix<TDomain>();
            std::string tag = bridge::GetDomainTag<TDomain>();
        }

        template<int dim>
        static void Dimension(Registry &reg, std::string grp) {
            std::string suffix = bridge::GetDimensionSuffix<dim>();
            std::string tag = bridge::GetDimensionTag<dim>();

        }

        template<typename TAlgebra>
        static void Algebra(Registry &reg, std::string grp) {
            std::string suffix = bridge::GetAlgebraSuffix<TAlgebra>();
            std::string tag = bridge::GetAlgebraTag<TAlgebra>();
        }

        // Memory Functions
        static void Common(Registry &reg, std::string grp) {
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
}

extern "C"
void InitUGPlugin_XBraidForUG4(Registry *reg, std::string param_grp) {
        using namespace xbraid;
        std::string grp = param_grp;

        grp.append("XBraidForUG4");

        // Space Time Communicator
        {
            using T_SpaceTimeCommunicator = SpaceTimeCommunicator ;
            const std::string name = "SpaceTimeCommunicator";
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
                    .add_method("set_openmp", &T_SpaceTimeCommunicator::set_openmp)
                    .set_construct_as_smart_pointer(true);
        }
        // ParallelLogger
        {
            using T_ParallelLogger = ParallelLogger ;
            const std::string name = "Paralog";
            reg->add_class_<T_ParallelLogger>(name, "XBraid", "")
                    .add_constructor()
                    .add_method("set_filename", &T_ParallelLogger::set_filename, "", "", "")
                    .add_method("set_comm", &T_ParallelLogger::set_comm, "", "", "")
                    .add_method("init", &T_ParallelLogger::init, "", "", "")
                    .add_method("release", &T_ParallelLogger::release, "", "", "")
                    .add_method("write", &T_ParallelLogger::write, "", "", "")
                    .set_construct_as_smart_pointer(true);
        }
        // ReplaceStandardStream
        {
            using T_ReplaceStandardStream = ReplaceStandardStream ;
            const std::string name = "ReplaceStandardStream";
            reg->add_class_<T_ReplaceStandardStream>(name, "XBraid", "")
                    .add_constructor()
                    .add_method("apply", &T_ReplaceStandardStream::apply, "", "", "")
                    .add_method("undo", &T_ReplaceStandardStream::undo, "", "", "")
                    .add_method("set_space_time_comm", &T_ReplaceStandardStream::set_space_time_comm, "", "", "")
                    .set_construct_as_smart_pointer(true);
        }
        // DiscardStandardStream
        {
        using T_DiscardStandardStream = DiscardStandardStream ;
        const std::string name = "DiscardStandardStream";
        reg->add_class_<T_DiscardStandardStream>(name, "XBraid", "")
                .add_constructor()
                .add_method("apply", &T_DiscardStandardStream::apply, "", "", "")
                .add_method("undo", &T_DiscardStandardStream::undo, "", "", "")
                .set_construct_as_smart_pointer(true);
       }
        // BraidTimer
        {
            using T_BraidTimer = BraidTimer ;
            reg->add_class_<T_BraidTimer>("BraidTimer", grp, "")
                    .add_constructor()
                    .add_method("start", &T_BraidTimer::start)
                    .add_method("stop", &T_BraidTimer::stop)
                    .add_method("get", &T_BraidTimer::get)
                    .set_construct_as_smart_pointer(true);
        }

    try {
        ug::bridge::RegisterCommon<Functionality>(*reg, grp);
        ug::bridge::RegisterDimensionDependent<Functionality>(*reg, grp);
        ug::bridge::RegisterDomainDependent<Functionality>(*reg, grp);
        ug::bridge::RegisterAlgebraDependent<Functionality>(*reg, grp);
        ug::bridge::RegisterDomainAlgebraDependent<Functionality>(*reg, grp);
    }
    UG_REGISTRY_CATCH_THROW(grp);

#ifdef XBraidPoroelasticity
        InitUGPlugin_XBraid_Poroelasticity(reg, grp);
#endif

}}