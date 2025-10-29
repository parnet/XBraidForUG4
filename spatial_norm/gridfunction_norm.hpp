#ifndef UGPLUGIN_XBRAIDFORUG4_GRIDFUNCTION_NORM_H
#define UGPLUGIN_XBRAIDFORUG4_GRIDFUNCTION_NORM_H

#include "interface/spatial_norm.hpp"



namespace ug{ namespace xbraid {

    template <typename TDomain, typename TAlgebra>
    class GridFunctionNorm : public BraidSpatialNorm<TDomain, TAlgebra> {
    public:

        //--------------------------------------------------------------------------------------------------------------

        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;

        using T_GridFunctionSpace = IGridFunctionSpace<T_GridFunction> ;
        using SP_GridFunctionSpace = SmartPtr<T_GridFunctionSpace> ;

        //--------------------------------------------------------------------------------------------------------------

        GridFunctionNorm() = default;
        ~GridFunctionNorm() override = default;

        //--------------------------------------------------------------------------------------------------------------

        void add_norm(SP_GridFunctionSpace norm){
            _norm = norm;
        }

        double norm(SP_GridFunction u) override {
            return _norm->norm2(*u);
        }

        //--------------------------------------------------------------------------------------------------------------

        SP_GridFunctionSpace _norm = SPNULL;
        //--------------------------------------------------------------------------------------------------------------
    };

}}
#endif