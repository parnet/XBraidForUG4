#ifndef UGPLUGIN_XBRAIDFORUG4_TRANSFER_RESTRICT_ELEMENTWISE_HPP
#define UGPLUGIN_XBRAIDFORUG4_TRANSFER_RESTRICT_ELEMENTWISE_HPP

#include "interface/restrict.hpp"



/* 2025-03-25

namespace ug{ namespace xbraid {

    template <typename TDomain, typename TAlgebra>
    class RestrictElementwiseWrapper : public IRestrict<TDomain,TAlgebra>{
    public:
        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;

        RestrictElementwiseWrapper() = default;

        ~RestrictElementwiseWrapper() override = default;

        void apply(SP_GridFunction cu, SP_GridFunction fu) override {

            const GridLevel& grid_level = fu->grid_level();

            T_GridFunction& uFine = *fu;
            T_GridFunction& uCoarse = *cu;
            RestrictElemwise(uCoarse,uFine);
        };

    };
}}*/

#endif