// todo generalize and split up functionalities
// interpolate result
// calculate norms
// write out in files
// more than one cmp?

#ifndef UGPLUGIN_XBRAIDFORUG4_OBSERVER_USERDATA_OBSERVER_HPP
#define UGPLUGIN_XBRAIDFORUG4_OBSERVER_USERDATA_OBSERVER_HPP


#include "lib_disc/spatial_disc/user_data/user_data.h"

#include "interface/observer_xbraid.hpp"
//2025-03 #include "../math/vector.h"

namespace ug{ namespace xbraid {

    template<typename TDomain, typename TAlgebra>
    class UserdataProcessEvaluateObserver : public IXBraidTimeIntegratorObserver<TDomain, TAlgebra> {
    public:
        //--------------------------------------------------------------------------------------------------------------

        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;

        using T_IDomainDiscretization = IDomainDiscretization<TAlgebra> ;
        using SP_DomainDisc = SmartPtr<T_IDomainDiscretization> ;

        using T_UserData = UserData<double, T_GridFunction::dim> ;
        using SP_UserData = SmartPtr<T_UserData> ;

        //--------------------------------------------------------------------------------------------------------------

        SP_UserData m_data;
        SP_DomainDisc m_domainDisc;

        std::ofstream *outfile{};

        const char *m_cmp{};
        bool m_relative = false;

        //--------------------------------------------------------------------------------------------------------------

        UserdataProcessEvaluateObserver() = default;

        ~UserdataProcessEvaluateObserver() override {
            outfile->close();
            delete outfile;
        }

        //--------------------------------------------------------------------------------------------------------------

        void set_domain(SP_DomainDisc p_domainDisc) {
            m_domainDisc = p_domainDisc;
        }

        void set_domain_disc(SP_DomainDisc p_domainDisc) { // todo delete one method
            m_domainDisc = p_domainDisc;
        }

        void set_relative(bool val) {
            this->m_relative = val;
        }

        void set_file(const char *filename) {
            outfile = new std::ofstream();
            outfile->open(filename);
            (*this->outfile) << "iteration; level; index; time; || u - v ||;";
            if (this->m_relative) {
                (*this->outfile) << "  || u ||; relative;";
            }
            (*this->outfile) << std::endl;
        }

        void set_generator_component(const char *cmp) {
            this->m_cmp = cmp;
        }

        void set_vector_generator(SP_UserData p_data) {
            this->m_data = p_data;
        }

        void set_vector_generator(const char *fctName) {
            set_vector_generator(LuaUserDataFactory<double, TDomain::dim>::create(fctName));
        }

        void set_vector_generator(LuaFunctionHandle fct) {
            set_vector_generator(make_sp(new LuaUserData<double, TDomain::dim>(fct)));
        }

        void write_time_pvd(SP_GridFunction u) {

        };

        bool step_process(SP_GridFunction u, int index, double time,double dt) override {
            return this->write(u, index, time, 0);
        };

        bool step_process(SP_GridFunction u, int index, double time, double dt,int iteration, int level) override {
            double relative;
            double unorm;

            // compute solution
            SP_GridFunction vec = u->clone_without_values();
            Interpolate(this->m_data, vec, this->m_cmp, NULL, time);
            m_domainDisc->adjust_solution(*vec.get(), time);

            // calculate difference / error
            ::ug::VecSubtract(*vec.get(), *vec.get(),  *u.get());

            // calculate absolut and relative norm
            double vecnorm = vec->norm();

            if (this->m_relative) {
                SP_GridFunction uvec = u->clone();
                unorm = uvec->norm();
                relative =  (vecnorm / unorm);
            }

            if (iteration == -1) {
                (*this->outfile) << std::setw(4) << "" << ";" << std::setw(4) << "" << ";";
            } else {
                (*this->outfile) << std::setw(4) << iteration << ";" << std::setw(4) << level << ";";
            }

            (*this->outfile) << std::setw(7) << index << ";"
                             << std::setw(12) << time << ";"
                             << std::setw(12) << vecnorm << ";";

            if (this->m_relative) {
                (*this->outfile) << std::setw(12) << unorm << ";" << std::setw(12) << (relative) << ";";
            }
            (*this->outfile) << std::endl;
            return true;
        };
        //--------------------------------------------------------------------------------------------------------------
    };
}}


#endif