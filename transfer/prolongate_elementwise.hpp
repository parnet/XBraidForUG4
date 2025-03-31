#ifndef UGPLUGIN_XBRAIDFORUG4_TRANSFER_ELEMENTWISE_HPP
#define UGPLUGIN_XBRAIDFORUG4_TRANSFER_ELEMENTWISE_HPP

#include "interface/prolongate.hpp"
#include "lib_disc/function_spaces/level_transfer.h"


/* 2025-03-25
namespace ug{ namespace xbraid {
    template <typename TDomain, typename TAlgebra>
    class ProlongateElementwiseWrapper: public IProlongate<TDomain,TAlgebra>  {
    public:
        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;

        ProlongateElementwiseWrapper() = default;

        ~ProlongateElementwiseWrapper() override = default;


        virtual void apply(SP_GridFunction fu, SP_GridFunction cu) {
            T_GridFunction& uFine = *fu;
            T_GridFunction& uCoarse = *cu;
            ProlongateElemwise(uFine, uCoarse);
        };
    };
}}
*/
#endif