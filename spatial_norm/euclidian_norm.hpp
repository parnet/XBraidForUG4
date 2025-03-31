#ifndef UGPLUGIN_XBRAIDFORUG4_SPATIAL_NORM_EUCLIDIAN_NORM_H
#define UGPLUGIN_XBRAIDFORUG4_SPATIAL_NORM_EUCLIDIAN_NORM_H

#include "interface/spatial_norm.hpp"





namespace ug{ namespace xbraid {

    template <typename TDomain, typename TAlgebra>
    class BraidEuclidianNorm : public BraidSpatialNorm<TDomain, TAlgebra> {
    public:

        //--------------------------------------------------------------------------------------------------------------

        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;

        //--------------------------------------------------------------------------------------------------------------

        BraidEuclidianNorm() = default;
        ~BraidEuclidianNorm() override = default;

        //--------------------------------------------------------------------------------------------------------------

        double norm(SP_GridFunction u) override {
            return u->norm();
        }

        //--------------------------------------------------------------------------------------------------------------
    };

}}
#endif