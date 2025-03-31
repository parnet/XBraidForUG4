#ifndef UGPLUGIN_XBRAIDFORUG4_POROELASTICITY_BIOT_BRAID_DISPLACEMENT_NORM_HPP
#define UGPLUGIN_XBRAIDFORUG4_POROELASTICITY_BIOT_BRAID_DISPLACEMENT_NORM_HPP


//2025-03 #include "lib_disc/function_spaces/integrate.h"
#include "lib_disc/function_spaces/grid_function.h"

#include "interface/spatial_norm.hpp"
#include "util/parallel_logger.hpp"

#include "biot_error_data.hpp"



namespace ug {  namespace xbraid { namespace poro {

        template<typename TDomain, typename TAlgebra>
        class BiotBraidDisplacementNorm : public BraidSpatialNorm<TDomain, TAlgebra> {
        public:

            //--------------------------------------------------------------------------------------------------------------

            using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
            using SP_GridFunction = SmartPtr<T_GridFunction> ;

            using T_ParallelLogger = ParallelLogger ;
            using SP_Paralog = SmartPtr<T_ParallelLogger> ;

            //--------------------------------------------------------------------------------------------------------------

                BiotBraidDisplacementNorm() : BraidSpatialNorm<TDomain, TAlgebra>() {}

            ~BiotBraidDisplacementNorm() override = default;

            //--------------------------------------------------------------------------------------------------------------

            int count = 0;

            SP_Paralog m_log;

            //--------------------------------------------------------------------------------------------------------------

            void set_log(SP_Paralog log) {
                this->m_log = log;
            }

            double norm(SP_GridFunction u) override {
                BiotErrorData<TDomain,TAlgebra> errdata = BiotErrorData<TDomain,TAlgebra>();
                errdata.compute(u);

                m_log->o << ">R> rnorm idx=" << count << std::endl;
                m_log->o << std::setw(10) << ">R>  l2(p)"
                         << std::setw(20) << errdata.l2_norm_p
                         << std::endl;

                m_log->o << std::setw(10) << ">R> l2(ux)"
                         << std::setw(20) << errdata.l2_norm_ux
                         << std::endl;

                m_log->o << std::setw(10) << ">R> l2(uy)"
                         << std::setw(20) << errdata.l2_norm_uy
                         << std::endl;

                m_log->o << std::setw(10) << ">R> h1(ux)"
                         << std::setw(20) << errdata.h1_norm_ux
                         << std::endl;

                m_log->o << std::setw(10) << ">R> h1(uy)"
                         << std::setw(20) << errdata.h1_norm_uy
                         << std::endl;

                count ++ ;
                return errdata.l2_norm_ux;
            }
            //--------------------------------------------------------------------------------------------------------------
        };
}}}
#endif