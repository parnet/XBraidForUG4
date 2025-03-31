#ifndef UGPLUGIN_XBRAIDFORUG4_POROELASTICITY_BRIDGE_HPP
#define UGPLUGIN_XBRAIDFORUG4_POROELASTICITY_BRIDGE_HPP
#include "bridge/util.h"
#include "bridge/util_domain_algebra_dependent.h"


#include "biot_error_data.hpp"
#include "braid_biot_control.hpp"
#include "braid_biot_estimator.hpp"
#include "braid_biot_precomputed.hpp"
#include "biot_braid_displacement_norm.hpp"
//2025-03 #include "braid_heat_check.hpp"

using namespace std;
using namespace ::ug::bridge;

namespace ug {
    namespace xbraid {
        namespace poro {
            struct Functionality {


                template<typename TDomain, typename TAlgebra>
                static void DomainAlgebra(Registry &reg, string grp) {

                    string suffix = GetDomainAlgebraSuffix<TDomain, TAlgebra>();
                    string tag = GetDomainAlgebraTag<TDomain, TAlgebra>();

                    { // BiotBraidDisplacementNorm
                        using T_BiotSpatialNorm = BiotBraidDisplacementNorm<TDomain, TAlgebra> ;
                        using T_SpatialNorm = BraidSpatialNorm<TDomain, TAlgebra> ;

                        string name = string("BiotBraidDisplacementNorm").append(suffix);
                        reg.add_class_<T_BiotSpatialNorm, T_SpatialNorm>(name, grp)
                                .add_constructor()
                                .add_method("norm", &T_BiotSpatialNorm::norm, "", "", "")
                                .add_method("set_log", &T_BiotSpatialNorm::set_log, "", "", "")
                                .set_construct_as_smart_pointer(true);
                        reg.add_class_to_group(name, "BiotBraidDisplacementNorm", tag);
                    }

                    { // BiotBraidDisplacementNorm
                        using T_BiotErrorData = BiotErrorData<TDomain, TAlgebra> ;

                        string name = string("BiotErrorData").append(suffix);
                        reg.add_class_<T_BiotErrorData>(name, grp)
                                .add_constructor()
                                .add_method("compute", &T_BiotErrorData::compute, "", "", "")
                                .add_method("set_order", &T_BiotErrorData::set_order, "", "", "")
                                .set_construct_as_smart_pointer(true);
                        reg.add_class_to_group(name, "BiotErrorData", tag);
                    }

                    {// BiotBraidSpatialNorm
                        using T_BiotSpatialNorm = BiotBraidSpatialNorm<TDomain, TAlgebra> ;
                        using T_SpatialNorm = BraidSpatialNorm<TDomain, TAlgebra> ;

                        string name = string("BiotBraidSpatialNorm").append(suffix);
                        reg.add_class_<T_BiotSpatialNorm, T_SpatialNorm>(name, grp)
                                .add_constructor()
                                .add_method("norm", &T_BiotSpatialNorm::norm, "", "", "")
                                .add_method("set_order", &T_BiotSpatialNorm::set_order, "", "", "")
                                .add_method("set_parameter", &T_BiotSpatialNorm::set_parameter, "", "", "")
                                .set_construct_as_smart_pointer(true);
                        reg.add_class_to_group(name, "BiotBraidSpatialNorm", tag);

                    }

                    {   //BraidBiotCheckPrecomputed
                        using T_BraidBiotCheckPrecomputed = BraidBiotCheckPrecomputed<TDomain, TAlgebra> ;
                        using T_IXBraidTimeIntegratorObserver = IXBraidTimeIntegratorObserver<TDomain, TAlgebra> ;
                        string name = string("BraidBiotCheckPrecomputed").append(suffix);
                        reg.add_class_<T_BraidBiotCheckPrecomputed, T_IXBraidTimeIntegratorObserver>(name, grp)
                                .add_constructor()
                                .add_method("compare_norms", &T_BraidBiotCheckPrecomputed::compare_norms, "", "", "")
                                .add_method("set_num_ref", &T_BraidBiotCheckPrecomputed::set_num_ref, "", "", "")
                                .add_method("set_max_index", &T_BraidBiotCheckPrecomputed::set_max_index, "", "", "")
                                .add_method("set_base_path", &T_BraidBiotCheckPrecomputed::set_base_path, "", "", "")
                                .add_method("set_c_factor", &T_BraidBiotCheckPrecomputed::set_c_factor, "", "", "")
                                .add_method("set_solution_name", &T_BraidBiotCheckPrecomputed::set_solution_name, "", "", "")
                                .add_method("set_diff_name", &T_BraidBiotCheckPrecomputed::set_diff_name, "", "", "")
                                .add_method("step_process", &T_BraidBiotCheckPrecomputed::lua_write, "", "", "")
                                .add_method("lua_compare", &T_BraidBiotCheckPrecomputed::lua_compare, "", "", "")
                                .add_method("set_log", &T_BraidBiotCheckPrecomputed::set_log, "", "", "")
                                .add_method("print", &T_BraidBiotCheckPrecomputed::print, "", "", "")
                                .add_method("set_vtk_write_mode", &T_BraidBiotCheckPrecomputed::set_vtk_write_mode, "", "", "")
                                .add_method("set_io_write_mode", &T_BraidBiotCheckPrecomputed::set_io_write_mode, "","", "")
                                .set_construct_as_smart_pointer(true);
                        reg.add_class_to_group(name, "BraidBiotCheckPrecomputed", tag);
                    }


                    /*{   //BraidHeatCheck
                        using T_BraidHeatCheck = BraidHeatCheck<TDomain, TAlgebra> ;
                        using T_IXBraidTimeIntegratorObserver = IXBraidTimeIntegratorObserver<TDomain, TAlgebra> ;
                        string name = string("BraidHeatCheck").append(suffix);
                        reg.add_class_<T_BraidHeatCheck, T_IXBraidTimeIntegratorObserver>(name, grp)
                                .add_constructor()
                                .add_method("set_num_ref", &T_BraidHeatCheck::set_num_ref, "", "", "")
                                .add_method("set_max_index", &T_BraidHeatCheck::set_max_index, "", "", "")
                                .add_method("set_c_factor", &T_BraidHeatCheck::set_c_factor, "", "", "")
                                .add_method("set_solution_name", &T_BraidHeatCheck::set_solution_name, "", "", "")
                                .add_method("set_diff_name", &T_BraidHeatCheck::set_diff_name, "", "", "")
                                .add_method("step_process", &T_BraidHeatCheck::lua_write, "", "", "")
                                .add_method("lua_compare", &T_BraidHeatCheck::lua_compare, "", "", "")
                                .add_method("set_log", &T_BraidHeatCheck::set_log, "", "", "")
                                .add_method("print", &T_BraidHeatCheck::print, "", "", "")
                                .add_method("set_vtk_write_mode", &T_BraidHeatCheck::set_vtk_write_mode, "", "", "")
                                .add_method("set_io_write_mode", &T_BraidHeatCheck::set_io_write_mode, "","", "")
                                .set_construct_as_smart_pointer(true);
                        reg.add_class_to_group(name, "BraidHeatCheck", tag);
                    }*/
                }


                template<typename TDomain>
                static void Domain(Registry &reg, string grp) {}

                template<int dim>
                static void Dimension(Registry &reg, string grp) {}

                template<typename TAlgebra>
                static void Algebra(Registry &reg, string grp) {}

                static void Common(Registry &reg, string grp) {}

            };

            struct FunctionalityFor2D {
                template<typename TDomain, typename TAlgebra>
                static void DomainAlgebra(Registry &reg, string grp) {
                    const string suffix = GetDomainAlgebraSuffix<TDomain, TAlgebra>();
                    const string tag = GetDomainAlgebraTag<TDomain, TAlgebra>();

                    // Braid Time Integrator
                    {
                        using T_BraidBiotCheck= BraidBiotCheck<TDomain, TAlgebra> ;
                        using T_SpatialNorm=  xbraid::IXBraidTimeIntegratorObserver<TDomain, TAlgebra> ;

                        using T_Problem = Poroelasticity::BarryMercerProblem<TDomain,TAlgebra> ;
                        using SP_Problem = SmartPtr<T_Problem> ;
                        string name = string("BraidBiotCheck").append(suffix);
                        reg.add_class_<T_BraidBiotCheck, T_SpatialNorm>(name, grp)
                                .add_constructor()
                                .add_method("set_problem", reinterpret_cast<void (T_BraidBiotCheck::*)(SP_Problem)>(&T_BraidBiotCheck::set_problem), "", "", "")
                                //.add_method("set_napprox", &T_BraidBiotCheck::set_napprox, "", "", "")
                                .add_method("step_process", &T_BraidBiotCheck::lua_write, "", "", "")
                                .add_method("set_filename", &T_BraidBiotCheck::set_filename, "", "", "")
                                .set_construct_as_smart_pointer(true);
                        reg.add_class_to_group(name, "BraidBiotCheck", tag);
                    }
                }
            };
        }



        void InitUGPlugin_XBraid_Poroelasticity(Registry *reg, string grp) {
            using namespace ::ug::xbraid::poro;
            grp.append("XBraidPoroelasticity");
            try {
                RegisterCommon<Functionality>(*reg, grp);
                RegisterDimensionDependent<Functionality>(*reg, grp);
                RegisterDomainDependent<Functionality>(*reg, grp);
                RegisterAlgebraDependent<Functionality>(*reg, grp);
                RegisterDomainAlgebraDependent<Functionality>(*reg, grp);
                RegisterDomain2dAlgebraDependent<FunctionalityFor2D>(*reg,grp);
            }
            UG_REGISTRY_CATCH_THROW(grp);

        }
    }
}
#endif