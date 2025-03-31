
#include "bridge/util.h"


#include "transfer/__package.hpp"

namespace ug{ namespace xbraid {
        struct RestrictProlongate_Functionality {

        template<typename TDomain, typename TAlgebra>
        static void DomainAlgebra(bridge::Registry &reg, std::string grp) {
            std::string suffix = bridge::GetDomainAlgebraSuffix<TDomain, TAlgebra>();
            std::string tag = bridge::GetDomainAlgebraTag<TDomain, TAlgebra>();
            /* 2025-03-25{
                using T_IProlongate = IProlongate<TDomain, TAlgebra> ;
                std::string name = std::string("IProlongate").append(suffix);
                reg.add_class_<T_IProlongate>(name,grp)
                        .add_method("apply", &T_IProlongate::apply, "", "", "")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name, "IProlongate", tag);
            }*//* 2025-03-25            {
                using T_RealCoarsen = RealCoarsen<TDomain, TAlgebra>;
                std::string name = std::string("RealCoarsen").append(suffix);
                reg.add_class_<T_RealCoarsen>(name,grp)
                        .add_constructor()
                        .add_method("set_transfer", &T_RealCoarsen::set_transfer, "", "", "")
                        .add_method("apply", &T_RealCoarsen::apply, "", "", "")
                        .add_method("make_nontop", &T_RealCoarsen::make_nontop, "", "", "")
                        .add_method("injection", &T_RealCoarsen::injection, "", "", "")
                        .add_method("prolong", &T_RealCoarsen::prolong, "", "", "")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name, "RealCoarsen", tag);
            } */

            /* 2025-03-25{
                using T_Prolongate = ProlongateWrapper<TDomain, TAlgebra> ;
                using T_IProlongate = IProlongate<TDomain, TAlgebra> ;
                std::string name = std::string("ProlongateWrapper").append(suffix);
                reg.add_class_<T_Prolongate,T_IProlongate>(name,grp)
                        .add_constructor()
                        .add_method("apply", &T_Prolongate::apply, "", "", "")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name, "ProlongateWrapper", tag);
            }*/
            /* 2025-03-25{
                using T_Prolongate = ProlongateElementwiseWrapper<TDomain, TAlgebra> ;
                using T_IProlongate = IProlongate<TDomain, TAlgebra> ;
                std::string name = std::string("ProlongateElementwiseWrapper").append(suffix);
                reg.add_class_<T_Prolongate,T_IProlongate>(name,grp)
                        .add_constructor()
                        .add_method("apply", &T_Prolongate::apply, "", "", "")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name, "ProlongateElementwiseWrapper", tag);
            }*/
            /* 2025-03-25{
                using T_Prolongate = ProlongateP1Wrapper<TDomain, TAlgebra> ;
                using T_IProlongate = IProlongate<TDomain, TAlgebra> ;
                std::string name = std::string("ProlongateP1Wrapper").append(suffix);
                reg.add_class_<T_Prolongate,T_IProlongate>(name,grp)
                        .add_constructor()
                        .add_method("apply", &T_Prolongate::apply, "", "", "")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name, "ProlongateP1Wrapper", tag);
            }*/
            /* 2025-03-25{
                using T_IRestrict = IRestrict<TDomain, TAlgebra> ;
                std::string name = std::string("IRestrict").append(suffix);
                reg.add_class_<T_IRestrict>(name,grp)
                        .add_method("apply", &T_IRestrict::apply, "", "", "")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name, "IRestrict", tag);
            }*/ /* 2025-03-25 {
                using T_Restrict = RestrictWrapper<TDomain, TAlgebra> ;
                using T_IRestrict = IRestrict<TDomain, TAlgebra> ;
                std::string name = std::string("RestrictWrapper").append(suffix);
                reg.add_class_<T_Restrict,T_IRestrict>(name,grp)
                        .add_constructor()
                        .add_method("apply", &T_Restrict::apply, "", "", "")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name, "RestrictWrapper", tag);
            } */
            /* 2025-03-25 {
                using T_Restrict = RestrictElementwiseWrapper<TDomain, TAlgebra> ;
                using T_IRestrict = IRestrict<TDomain, TAlgebra> ;
                std::string name = std::string("RestrictElementwiseWrapper").append(suffix);
                reg.add_class_<T_Restrict,T_IRestrict>(name,grp)
                        .add_constructor()
                        .add_method("apply", &T_Restrict::apply, "", "", "")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name, "RestrictElementwiseWrapper", tag);
            }*/ /* 2025-03-25   {
                using T_Restrict = RestrictP1Wrapper<TDomain, TAlgebra> ;
                using T_IRestrict = IRestrict<TDomain, TAlgebra> ;
                std::string name = std::string("RestrictP1Wrapper").append(suffix);
                reg.add_class_<T_Restrict,T_IRestrict>(name,grp)
                        .add_constructor()
                        .add_method("apply", &T_Restrict::apply, "", "", "")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name, "RestrictP1Wrapper", tag);
            }*/

            }
        };
    }
}

