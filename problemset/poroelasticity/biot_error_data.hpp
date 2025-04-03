#ifndef UGPLUGIN_XBRAIDFORUG4_POROELASTICITY_BIOT_ERROR_DATA_HPP
#define UGPLUGIN_XBRAIDFORUG4_POROELASTICITY_BIOT_ERROR_DATA_HPP

//2025-03 #include <math.h>

#include "lib_disc/function_spaces/integrate.h"
#include "lib_disc/function_spaces/grid_function.h"

namespace ug { namespace xbraid { namespace poro {

    template<typename TDomain, typename TAlgebra>
    class BiotErrorData {
    public:

        //--------------------------------------------------------------------------------------------------------------

        using T_GridFunction = GridFunction<TDomain, TAlgebra>;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;


        //--------------------------------------------------------------------------------------------------------------

        void compute(SP_GridFunction u) {

            this->l2_norm_p_ = L2Norm(*u.get(), "p", porder_);;
            this->l2_norm_ux_ = L2Norm(*u.get(), "ux", uorder_);
            this->l2_norm_uy_ = L2Norm(*u.get(), "uy", uorder_);

            this->h1_norm_ux_ = H1SemiNorm(*u.get(), "ux", uorder_);
            this->h1_norm_uy_ = H1SemiNorm(*u.get(), "uy", uorder_);

        }

        void set_order(int porder = 2, int uorder = 4) {
            this->porder_ = porder;
            this->uorder_ = uorder;
        }
        //--------------------------------------------------------------------------------------------------------------

        double l2_norm_p_ = 0.0;
        double l2_norm_ux_ = 0.0;
        double l2_norm_uy_ = 0.0;

        double h1_norm_ux_ = 0.0;
        double h1_norm_uy_ = 0.0;

        int porder_ = 2;
        int uorder_ = 4;

        //--------------------------------------------------------------------------------------------------------------
    };

}}}

#endif