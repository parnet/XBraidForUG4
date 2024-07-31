#ifndef UG_PLUGIN_XBRAIDFORUG4_ZERO_INITIALIZER_H
#define UG_PLUGIN_XBRAIDFORUG4_ZERO_INITIALIZER_H

#include "../interface/initializer.h"


namespace ug { namespace XBraidForUG4 {
    template <typename TDomain, typename TAlgebra>
    class ZeroInitializer : public BraidInitializer<TDomain, TAlgebra> {
    public:
        //--------------------------------------------------------------------------------------------------------------

        typedef GridFunction<TDomain, TAlgebra> T_GridFunction;
        typedef SmartPtr<T_GridFunction> SP_GridFunction;

        //--------------------------------------------------------------------------------------------------------------

        ZeroInitializer() = default;

        ~ZeroInitializer() override = default;

        //--------------------------------------------------------------------------------------------------------------

        void initialize(SP_GridFunction& u, number time) override {
            if (time == this->t_start) {
                u = this->m_u0->clone();
            } else {
                u = this->m_u0->clone_without_values();
            }
        }

        //--------------------------------------------------------------------------------------------------------------
    };
}}
#endif //UG_PLUGIN_XBRAIDFORUG4_ZERO_INITIALIZER_H
