#ifndef UGPLUGIN_XBRAIDFORUG4_INTERFACE_SPATIAL_NORM_H
#define UGPLUGIN_XBRAIDFORUG4_INTERFACE_SPATIAL_NORM_H

#include "lib_disc/function_spaces/grid_function.h"





namespace ug{ namespace xbraid {

    template <typename TDomain, typename TAlgebra>
    class BraidSpatialNorm {
    public:

        //--------------------------------------------------------------------------------------------------------------

        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;

        //--------------------------------------------------------------------------------------------------------------
    protected:
        BraidSpatialNorm() = default;

    public:
        virtual ~BraidSpatialNorm() = default;

        //--------------------------------------------------------------------------------------------------------------

        virtual double norm(SP_GridFunction u) = 0;

        //--------------------------------------------------------------------------------------------------------------
    };


}}
#endif