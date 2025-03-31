#ifndef UGPLUGIN_XBRAIDFORUG4_POROELASTICITY_BRAID_BIOT_ESTIMATOR_HPP
#define UGPLUGIN_XBRAIDFORUG4_POROELASTICITY_BRAID_BIOT_ESTIMATOR_HPP


#include "lib_disc/function_spaces/integrate.h"
#include "lib_disc/function_spaces/grid_function.h"

#include "interface/spatial_norm.hpp"



namespace ug { namespace xbraid {namespace poro {

    template<typename TDomain, typename TAlgebra>
    class BiotBraidSpatialNorm final : public BraidSpatialNorm<TDomain, TAlgebra> {
    public:

        //--------------------------------------------------------------------------------------------------------------

        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;

        //--------------------------------------------------------------------------------------------------------------

        int m_uorder = 4;
        int m_porder = 2;
        double p_factor = 1;
        double u_factor = 1;


        //--------------------------------------------------------------------------------------------------------------

        BiotBraidSpatialNorm() : BraidSpatialNorm<TDomain, TAlgebra>() {}

        ~BiotBraidSpatialNorm() override = default;

        //--------------------------------------------------------------------------------------------------------------

        void set_order(int uorder, int porder) {
            this->m_uorder = uorder;
            this->m_porder = porder;
        }

        void set_parameter(double alpha, double lambda, double mu) {
            p_factor = alpha;
            u_factor = lambda + 2 * mu;
        }

        double norm(SP_GridFunction u) override {
            const double norm_x = H1SemiNorm(*u.get(), "ux", this->m_uorder);
            const double norm_y = H1SemiNorm(*u.get(), "uy", this->m_uorder);
            const double norm_p = L2Norm(*u.get(), "p", this->m_porder);


            const double pnorm = p_factor * norm_p * norm_p; //p_factor*unorm_p*unorm_p;
            const double unorm = u_factor * (norm_x * norm_x + norm_y * norm_y);

            const double total_norm = sqrt(unorm + pnorm);

            return total_norm;
        }
        //--------------------------------------------------------------------------------------------------------------
    };
}}}

#endif