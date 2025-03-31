#ifndef UGPLUGIN_XBRAIDFORUG4_TRANSFER_REAL_COARSEN_HPP
#define UGPLUGIN_XBRAIDFORUG4_TRANSFER_REAL_COARSEN_HPP
//2025-03 #include "spatial_grid_transfer.hpp"
//2025-03 #include "lib_disc/operator/linear_operator/std_transfer.h"

/* 2025-03-25
#include "util/pragma.hpp"
#include "lib_disc/operator/linear_operator/std_injection.h"



namespace ug { namespace xbraid {

    template<typename TDomain, typename TAlgebra>
    class RealCoarsen {
    public:
        using T_GridFunction = GridFunction<TDomain, TAlgebra>;
        using SP_GridFunction = SmartPtr<T_GridFunction>;

        SmartPtr<ITransferOperator<TDomain, TAlgebra>> transfer;

        int num_ref = 0;

        void set_transfer(SmartPtr<ITransferOperator<TDomain, TAlgebra> > transfer) {
            this->transfer = transfer;
        }


        SP_GridFunction make_nontop(SP_GridFunction sp_gf) {
            __debug(std::cout << "GridFunction was Top" << std::endl);
            SmartPtr<ApproximationSpace<TDomain> > approx_space = sp_gf->approx_space();
            const size_t level = approx_space->num_levels();
            __debug(std::cout << "Transforming to level = " << level << std::endl);
            SP_GridFunction grid_function = make_sp(new T_GridFunction(approx_space, level - 1, false));

            const size_t usize = sp_gf->size();
            for (size_t i = 0; i < usize; ++i) {
                (*grid_function)[i] = (*sp_gf)[i];
            }

            grid_function->enable_redistribution(sp_gf->redistribution_enabled());
            grid_function->set_storage_type(sp_gf->get_storage_mask());
            grid_function->set_layouts(sp_gf->layouts());

            return grid_function;
        }





        SP_GridFunction apply(SP_GridFunction sp_fine) {
            // transforms from top vector to non top if neccessary
            // create coarse target vector
            // do restrict
            // sets values from original vector

            SP_GridFunction sp_fine_nontop;
            // ---------------------------------------------------------------------------------------------------------
            //(std::cout << "surface: " << sp_fine->grid_level().is_surface() << std::endl << std::flush);
            //(std::cout << "ghosts: " << sp_fine->grid_level().ghosts() << std::endl << std::flush);
            //(std::cout << "level: " << sp_fine->grid_level().is_level() << std::endl << std::flush);
            //(std::cout << "level: " << sp_fine->grid_level().level() << std::endl << std::flush);
            //(std::cout << "top: " << sp_fine->grid_level().top() << std::endl << std::flush);
            //(std::cout << "type: " << sp_fine->grid_level().type() << std::endl << std::flush);
            // ---------------------------------------------------------------------------------------------------------
            size_t target_level;
            SmartPtr<ApproximationSpace<TDomain> > approx_space = sp_fine->approx_space();
            if (sp_fine->grid_level().level() == GridLevel::TOP) {
                sp_fine_nontop = this->make_nontop(sp_fine);
                target_level = sp_fine->approx_space()->num_levels() - 2;
            } else {
                sp_fine_nontop = sp_fine;
                target_level = sp_fine->grid_level().level() - 1;
            }
            // ---------------------------------------------------------------------------------------------------------
            //(std::cout << "surface: " << sp_fine_nontop->grid_level().is_surface() << std::endl << std::flush);
            //(std::cout << "ghosts: " << sp_fine_nontop->grid_level().ghosts() << std::endl << std::flush);
            //(std::cout << "level: " << sp_fine_nontop->grid_level().is_level() << std::endl << std::flush);
            //(std::cout << "level: " << sp_fine_nontop->grid_level().level() << std::endl << std::flush);
            //(std::cout << "top: " << sp_fine_nontop->grid_level().top() << std::endl << std::flush);
            //(std::cout << "type: " << sp_fine_nontop->grid_level().type() << std::endl << std::flush);
            // ---------------------------------------------------------------------------------------------------------
            __debug(std::cout << "NumRef: " << target_level +1<< " ---> " << target_level<< " | "<< std::endl <<std::flush);

            T_GridFunction *gf_fu = sp_fine_nontop.get();
            auto *gf_cu = new T_GridFunction(approx_space, target_level, false);
            //this->restrict->apply(*sp_cu, sp_fu);
            this->transfer->do_restrict(*gf_cu, *gf_fu);

            SP_GridFunction sp_cu = make_sp(gf_cu);
            sp_cu->set_storage_type(sp_fine_nontop->get_storage_type());
            sp_cu->enable_redistribution(sp_fine_nontop->redistribution_enabled());
            //grid_function->set_storage_type(u->get_storage_mask());
            //sp_cu->set_layouts(u->layouts());
            return sp_cu;
        }

        SP_GridFunction injection(SP_GridFunction sp_fine) {
            // ---------------------------------------------------------------------------------------------------------
            //(std::cout << "surface: " << sp_fine->grid_level().is_surface() << std::endl << std::flush);
            //(std::cout << "ghosts: " << sp_fine->grid_level().ghosts() << std::endl << std::flush);
            //(std::cout << "level: " << sp_fine->grid_level().is_level() << std::endl << std::flush);
            //(std::cout << "level: " << sp_fine->grid_level().level() << std::endl << std::flush);
            //(std::cout << "top: " << sp_fine->grid_level().top() << std::endl << std::flush);
            //(std::cout << "type: " << sp_fine->grid_level().type() << std::endl << std::flush);
            // ---------------------------------------------------------------------------------------------------------
            SP_GridFunction sp_fine_nontop;
            size_t target_level;
            SmartPtr<ApproximationSpace<TDomain> > approx_space = sp_fine->approx_space();
            SmartPtr<StdInjection<TDomain,TAlgebra>> injection_operator = make_sp(new StdInjection<TDomain,TAlgebra>(approx_space));
            if (sp_fine->grid_level().level() == GridLevel::TOP) {
                sp_fine_nontop = this->make_nontop(sp_fine);
            } else {
                sp_fine_nontop = sp_fine;
            }

            // ---------------------------------------------------------------------------------------------------------
            (std::cout << "surface: " << sp_fine_nontop->grid_level().is_surface() << std::endl << std::flush);
            (std::cout << "ghosts: " << sp_fine_nontop->grid_level().ghosts() << std::endl << std::flush);
            (std::cout << "level: " << sp_fine_nontop->grid_level().is_level() << std::endl << std::flush);
            (std::cout << "level: " << sp_fine_nontop->grid_level().level() << std::endl << std::flush);
            (std::cout << "top: " << sp_fine_nontop->grid_level().top() << std::endl << std::flush);
            (std::cout << "type: " << sp_fine_nontop->grid_level().type() << std::endl << std::flush);
            // ---------------------------------------------------------------------------------------------------------
            target_level = sp_fine_nontop->grid_level().level() - 1;
            std::cout << "NumRef: " << target_level +1<< " ---> " << target_level<< " | "<< std::endl <<std::flush;
            auto *gf_cu = new T_GridFunction(approx_space, target_level, false);
            SP_GridFunction sp_cu = make_sp(gf_cu);

            std::cout << "NumRef: " << sp_fine_nontop->grid_level().level()<< " ---> " << sp_cu->grid_level().level() << " | "<< std::endl <<std::flush;

            injection_operator->set_levels(sp_cu->grid_level(), sp_fine_nontop->grid_level());
            injection_operator->init();
            T_GridFunction *gf_fu = sp_fine_nontop.get();
            injection_operator->do_restrict(*sp_cu, *sp_fine_nontop);
            sp_cu->set_storage_type(sp_fine_nontop->get_storage_type());
            sp_cu->enable_redistribution(sp_fine_nontop->redistribution_enabled());
            //grid_function->set_storage_type(u->get_storage_mask());
            //sp_cu->set_layouts(u->layouts());
            return sp_cu;
        }

        SP_GridFunction inject_prolong(SP_GridFunction sp_coarse) {
            // ---------------------------------------------------------------------------------------------------------
            //(std::cout << "surface: " << sp_coarse->grid_level().is_surface() << std::endl << std::flush);
            //(std::cout << "ghosts: " << sp_coarse->grid_level().ghosts() << std::endl << std::flush);
            //(std::cout << "level: " << sp_coarse->grid_level().is_level() << std::endl << std::flush);
            //(std::cout << "level: " << sp_coarse->grid_level().level() << std::endl << std::flush);
            //(std::cout << "top: " << sp_coarse->grid_level().top() << std::endl << std::flush);
            //(std::cout << "type: " << sp_coarse->grid_level().type() << std::endl << std::flush);
            // ---------------------------------------------------------------------------------------------------------
            size_t target_level;
            SmartPtr<ApproximationSpace<TDomain> > approx_space = sp_coarse->approx_space();
            SmartPtr<StdInjection<TDomain,TAlgebra>> injection_operator = make_sp(new StdInjection<TDomain,TAlgebra>(approx_space));


            // ---------------------------------------------------------------------------------------------------------
            target_level = sp_coarse->grid_level().level() + 1;;
            std::cout << "NumRef: " << target_level << " ---> " << target_level +1<< " | "<< std::endl <<std::flush;
            auto *gf_fine = new T_GridFunction(approx_space, target_level, false);
            SP_GridFunction sp_fine = make_sp(gf_fine);

            std::cout << "NumRef: " << sp_coarse->grid_level().level()<< " ---> " << sp_fine->grid_level().level() << " | "<< std::endl <<std::flush;

            injection_operator->set_levels(sp_coarse->grid_level(), sp_fine->grid_level());
            injection_operator->init();
            T_GridFunction *gf_fu = sp_fine.get();
            injection_operator->prolongate(*sp_fine, *sp_coarse);
            sp_fine->set_storage_type(sp_coarse->get_storage_type());
            sp_fine->enable_redistribution(sp_coarse->redistribution_enabled());
            //grid_function->set_storage_type(u->get_storage_mask());
            //sp_cu->set_layouts(u->layouts());
            return sp_fine;
        }

        SP_GridFunction prolong(SP_GridFunction sp_coarse) {
            // ---------------------------------------------------------------------------------------------------------
            //(std::cout << "surface: " << sp_coarse->grid_level().is_surface() << std::endl << std::flush);
            //(std::cout << "ghosts: " << sp_coarse->grid_level().ghosts() << std::endl << std::flush);
            //(std::cout << "level: " << sp_coarse->grid_level().is_level() << std::endl << std::flush);
            //(std::cout << "level: " << sp_coarse->grid_level().level() << std::endl << std::flush);
            //(std::cout << "top: " << sp_coarse->grid_level().top() << std::endl << std::flush);
            //(std::cout << "type: " << sp_coarse->grid_level().type() << std::endl << std::flush);
            // ---------------------------------------------------------------------------------------------------------
            size_t target_level;
            SmartPtr<ApproximationSpace<TDomain> > approx_space = sp_coarse->approx_space();
            SmartPtr<StdInjection<TDomain,TAlgebra>> injection_operator = make_sp(new StdInjection<TDomain,TAlgebra>(approx_space));


            // ---------------------------------------------------------------------------------------------------------
            target_level = sp_coarse->grid_level().level() + 1;;
            std::cout << "NumRef: " << target_level << " ---> " << target_level +1<< " | "<< std::endl <<std::flush;
            auto *gf_fine = new T_GridFunction(approx_space, target_level, false);
            SP_GridFunction sp_fine = make_sp(gf_fine);

            std::cout << "NumRef: " << sp_coarse->grid_level().level()<< " ---> " << sp_fine->grid_level().level() << " | "<< std::endl <<std::flush;

            this->transfer->set_levels(sp_coarse->grid_level(), sp_fine->grid_level());
            this->transfer->init();
            this->transfer->prolongate(*sp_fine, *sp_coarse);
            sp_fine->set_storage_type(sp_coarse->get_storage_type());
            sp_fine->enable_redistribution(sp_coarse->redistribution_enabled());
            //grid_function->set_storage_type(u->get_storage_mask());
            //sp_cu->set_layouts(u->layouts());
            return sp_fine;
        }
    };
}}
*/
#endif