#ifndef UGPLUGIN_XBRAIDFORUG4_TRANSFER_RESTRICT_P1_HPP
#define UGPLUGIN_XBRAIDFORUG4_TRANSFER_RESTRICT_P1_HPP

#include "interface/restrict.hpp"
#include "lib_disc/function_spaces/level_transfer.h"


/* 2025-03-25

namespace ug{ namespace xbraid {

    template <typename TDomain, typename TAlgebra>
    class RestrictP1Wrapper : public IRestrict<TDomain,TAlgebra> {
    public:
        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;

        RestrictP1Wrapper() = default;

        ~RestrictP1Wrapper() override = default;

        void apply(SP_GridFunction cu, SP_GridFunction fu) override{
            T_GridFunction& uFine = *fu;
            T_GridFunction& uCoarse = *cu;
            RestrictP1(uCoarse,uFine);
        };

    };
}}
*/
#endif