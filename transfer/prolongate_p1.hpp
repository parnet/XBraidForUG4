#ifndef UGPLUGIN_XBRAIDFORUG4_TRANSFER_PROLONGATE_P1_HPP
#define UGPLUGIN_XBRAIDFORUG4_TRANSFER_PROLONGATE_P1_HPP

#include "interface/prolongate.hpp"
#include "lib_disc/function_spaces/level_transfer.h"


/* 2025-03-25

namespace ug{ namespace xbraid {

template <typename TDomain, typename TAlgebra>
class ProlongateP1Wrapper : public IProlongate<TDomain,TAlgebra> {

public:
    using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
    using SP_GridFunction = SmartPtr<T_GridFunction> ;

    ProlongateP1Wrapper() = default;

    ~ProlongateP1Wrapper() override = default;


    void apply(SP_GridFunction fu, SP_GridFunction cu) override {
        T_GridFunction& uFine = *fu;
        T_GridFunction& uCoarse = *cu;
        ProlongateP1(uFine, uCoarse);
    };
};
}}
*/
#endif