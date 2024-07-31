// todo use a more general class

#ifndef UG_PLUGIN_XBRAIDFORUG4_SPATIAL_NORM_H
#define UG_PLUGIN_XBRAIDFORUG4_SPATIAL_NORM_H

#include "lib_disc/function_spaces/grid_function.h"


namespace ug { namespace XBraidForUG4 {
    template <typename TDomain, typename TAlgebra>
    class BraidSpatialNorm {
    public:

        //--------------------------------------------------------------------------------------------------------------

        typedef GridFunction<TDomain, TAlgebra> T_GridFunction;
        typedef SmartPtr<T_GridFunction> SP_GridFunction;

        //--------------------------------------------------------------------------------------------------------------

        BraidSpatialNorm() = default;

        virtual ~BraidSpatialNorm() = default;

        //--------------------------------------------------------------------------------------------------------------

        virtual double norm(SP_GridFunction u) = 0;

        //--------------------------------------------------------------------------------------------------------------
    };


}}
#endif //UG_PLUGIN_XBRAIDFORUG4_SPATIAL_NORM_H
