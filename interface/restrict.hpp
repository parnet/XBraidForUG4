#ifndef UGPLUGIN_XBRAIDFORUG4_TRANSFER_RESTRICT_INTERFACE_HPP
#define UGPLUGIN_XBRAIDFORUG4_TRANSFER_RESTRICT_INTERFACE_HPP

#include "common/util/smart_pointer.h"
#include "lib_disc/function_spaces/grid_function.h"



/* 2025-03-25
namespace ug{ namespace xbraid {

    template <typename TDomain, typename TAlgebra>
    class IRestrict {
    public:
        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;

    protected:
        IRestrict() = default;
    public:
        virtual ~IRestrict() = default;

        virtual void apply(SP_GridFunction cu, SP_GridFunction fu) = 0;


    };
}}*/

#endif