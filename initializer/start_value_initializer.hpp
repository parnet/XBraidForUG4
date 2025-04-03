#ifndef UGPLUGIN_XBRAIDFORUG4_INITIALIZER_START_VALUE_INITIALIZER_HPP
#define UGPLUGIN_XBRAIDFORUG4_INITIALIZER_START_VALUE_INITIALIZER_HPP

#include "interface/initializer.hpp"





namespace ug{ namespace xbraid {

    template <typename TDomain, typename TAlgebra>
    class GridFunctionInitializer : public BraidInitializer<TDomain, TAlgebra> {
    public:
        //--------------------------------------------------------------------------------------------------------------

        using T_GridFunction =  GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction>;

        //--------------------------------------------------------------------------------------------------------------

        GridFunctionInitializer() = default;

        ~GridFunctionInitializer() override = default;

        //--------------------------------------------------------------------------------------------------------------

        void initialize(SP_GridFunction& u, double time) override {
            u = this->u0_->clone();
        }

        //--------------------------------------------------------------------------------------------------------------
    };
}}

#endif