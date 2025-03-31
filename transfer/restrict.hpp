#ifndef UGPLUGIN_XBRAIDFORUG4_TRANSFER_RESTRICT_HPP
#define UGPLUGIN_XBRAIDFORUG4_TRANSFER_RESTRICT_HPP

#include "interface/restrict.hpp"
#include "lib_disc/function_spaces/level_transfer.h"



/* 2025-03-25
namespace ug{ namespace xbraid {

    template <typename TDomain, typename TAlgebra>
    class RestrictWrapper : public IRestrict<TDomain,TAlgebra>{
    public:
        using T_GridFunction=  GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;

        RestrictWrapper() = default;

        ~RestrictWrapper() override = default;

        void apply(SP_GridFunction cu, SP_GridFunction fu) override {
            T_GridFunction& uFine = *fu;
            T_GridFunction& uCoarse = *cu;

            std::cout << "Coarse = " << uCoarse.dof_distribution()->grid_level().level() << std::endl;
            std::cout << "  Fine = " << uFine.dof_distribution()->grid_level().level() << std::endl;

            Restrict(uCoarse,uFine);
        };

    };
}}*/

#endif