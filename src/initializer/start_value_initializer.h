#ifndef UG_PLUGIN_XBRAIDFORUG4_START_VALUE_INITIALIZER_H
#define UG_PLUGIN_XBRAIDFORUG4_START_VALUE_INITIALIZER_H

#include "../interface/initializer.h"


namespace ug { namespace XBraidForUG4 {
    template <typename TDomain, typename TAlgebra>
    class GridFunctionInitializer : public BraidInitializer<TDomain, TAlgebra> {
    public:
        //--------------------------------------------------------------------------------------------------------------

        typedef GridFunction<TDomain, TAlgebra> T_GridFunction;
        typedef SmartPtr<T_GridFunction> SP_GridFunction;

        //--------------------------------------------------------------------------------------------------------------

        GridFunctionInitializer() = default;

        ~GridFunctionInitializer() override = default;

        //--------------------------------------------------------------------------------------------------------------

        void initialize(SP_GridFunction& u, double time) override {
            u = this->m_u0->clone();
        }

        //--------------------------------------------------------------------------------------------------------------
    };
}}

#endif //UG_PLUGIN_XBRAIDFORUG4_START_VALUE_INITIALIZER_H
