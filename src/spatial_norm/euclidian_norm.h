//todo check
#ifndef UG_PLUGIN_XBRAIDFORUG4_EUCLIDIAN_NORM_H
#define UG_PLUGIN_XBRAIDFORUG4_EUCLIDIAN_NORM_H

#include "../interface/spatial_norm.h"



namespace ug { namespace XBraidForUG4 {

    template <typename TDomain, typename TAlgebra>
    class BraidEuclidianNorm : public BraidSpatialNorm<TDomain, TAlgebra> {
    public:

        //--------------------------------------------------------------------------------------------------------------

        typedef GridFunction<TDomain, TAlgebra> T_GridFunction;
        typedef SmartPtr<T_GridFunction> SP_GridFunction;

        //--------------------------------------------------------------------------------------------------------------

        BraidEuclidianNorm() : BraidSpatialNorm<TDomain, TAlgebra>() {}
        ~BraidEuclidianNorm() override = default;

        //--------------------------------------------------------------------------------------------------------------

        double norm(SP_GridFunction u) override {
            return u->norm();
        }

        //--------------------------------------------------------------------------------------------------------------
    };

}}
#endif //UG_PLUGIN_XBRAIDFORUG4_EUCLIDIAN_NORM_H
