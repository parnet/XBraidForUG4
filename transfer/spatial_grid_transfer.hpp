#ifndef UGPLUGIN_XBRAIDFORUG4_TRANSFER_SPATIALGRIDTRANSFER_HPP
#define UGPLUGIN_XBRAIDFORUG4_TRANSFER_SPATIALGRIDTRANSFER_HPP

#include <lib_disc/operator/linear_operator/std_injection.h>
#include <lib_disc/operator/linear_operator/std_transfer.h>
#include <lib_disc/operator/linear_operator/transfer_interface.h>
#include "lib_disc/spatial_disc/domain_disc.h"
#include "lib_disc/function_spaces/grid_function.h"

namespace ug{ namespace xbraid {

template <typename TDomain, typename TAlgebra>
class SpatialGridTransfer {
public:
    //--------------------------------------------------------------------------------------------------------------
    using T_GridFunction = GridFunction<TDomain,TAlgebra>;
    using T_ApproximationSpace = ApproximationSpace<TDomain>;
    using T_TransferOperator = ITransferOperator<TDomain,TAlgebra>;
    using T_DomainDisc = DomainDiscretization<TDomain,TAlgebra>;

    using SP_GridFunction = SmartPtr<T_GridFunction>;
    using SP_ApproximationSpace = SmartPtr<T_ApproximationSpace>;
    using SP_TransferOperator = SmartPtr<T_TransferOperator>;
    using SP_DomainDisc = SmartPtr<T_DomainDisc>;
    //--------------------------------------------------------------------------------------------------------------
    SpatialGridTransfer() = default;
    ~SpatialGridTransfer() = default;
    //--------------------------------------------------------------------------------------------------------------
    SP_GridFunction make_nontop(SP_GridFunction sp_gf);
    SP_GridFunction prolongate(SP_GridFunction sp_gf_coarse);
    SP_GridFunction restrict(SP_GridFunction sp_gf_fine);

    void set_approx_space(SP_ApproximationSpace approx_space);
    void set_prolongation(SP_TransferOperator prolongation);
    void set_restriction(SP_TransferOperator restriction);
    void set_domain(SP_DomainDisc domain_disc);
    void set_transfer(SP_TransferOperator transfer);
    void init();
    //--------------------------------------------------------------------------------------------------------------
protected:
private:
    //--------------------------------------------------------------------------------------------------------------
    SP_ApproximationSpace approximation_space_;
    SP_DomainDisc domain_disc_;

    SP_TransferOperator restriction_;
    //SP_TransferOperator andra_restriction_;
    SP_TransferOperator prolongation_;


    std::vector<SP_TransferOperator> level_restriction_;
    //std::vector<SP_TransferOperator> andra_level_restriction_;
    std::vector<SP_TransferOperator> level_prolongation_;

};


template<typename TDomain, typename TAlgebra>
SmartPtr<GridFunction<TDomain,TAlgebra>> SpatialGridTransfer<TDomain, TAlgebra>::make_nontop(SP_GridFunction sp_gf) {
    //__debug(std::cout << "GridFunction was Top" << std::endl);
    const size_t level = this->approximation_space_->num_levels();
    //__debug(std::cout << "Transforming to level = " << level - 1 << std::endl);
    SP_GridFunction grid_function = make_sp(new T_GridFunction(this->approximation_space_, level - 1, false));

    const size_t usize = sp_gf->size();
    for (size_t i = 0; i < usize; ++i) {
        (*grid_function)[i] = (*sp_gf)[i];
    }

    grid_function->enable_redistribution(sp_gf->redistribution_enabled());
    grid_function->set_storage_type(sp_gf->get_storage_mask());
    grid_function->set_layouts(sp_gf->layouts());

    return grid_function;
}

template<typename TDomain, typename TAlgebra>
SmartPtr<GridFunction<TDomain,TAlgebra>>
    SpatialGridTransfer<TDomain, TAlgebra>::prolongate(SP_GridFunction sp_gf_coarse) {
    // ---------------------------------------------------------------------------------------------------------
            size_t target_level = sp_gf_coarse->grid_level().level() + 1;;
            auto *gf_fine = new T_GridFunction(this->approximation_space_, target_level, false);
            SP_GridFunction sp_fine = make_sp(gf_fine);
            //std::cout << "ø apply prolongate to target level : " << target_level << std::endl;
            this->level_prolongation_[target_level]->prolongate(*sp_fine, *sp_gf_coarse);

            sp_fine->set_storage_type(sp_gf_coarse->get_storage_type());
            sp_fine->enable_redistribution(sp_gf_coarse->redistribution_enabled());

            return sp_fine;
}

template<typename TDomain, typename TAlgebra>
SmartPtr<GridFunction<TDomain,TAlgebra>> SpatialGridTransfer<TDomain, TAlgebra>::restrict(SP_GridFunction sp_gf_fine) {
            SP_GridFunction sp_fine_nontop;
            SmartPtr<ApproximationSpace<TDomain> > approx_space = sp_gf_fine->approx_space();
            if (sp_gf_fine->grid_level().level() == GridLevel::TOP) {
                sp_fine_nontop = this->make_nontop(sp_gf_fine);
            } else {
                sp_fine_nontop = sp_gf_fine;
            }

            size_t target_level = sp_fine_nontop->grid_level().level() - 1;
            //std::cout << "ø apply restrict to target level : " << target_level << std::endl;
            auto *gf_cu = new T_GridFunction(approx_space, target_level, false);
            SP_GridFunction sp_cu = make_sp(gf_cu);
            this->level_restriction_[target_level+1]->do_restrict(*sp_cu, *sp_fine_nontop);

            //auto *andra_gf_cu = new T_GridFunction(approx_space, target_level, false);
            //SP_GridFunction sp_andra_cu = make_sp(andra_gf_cu);
            //this->andra_level_restriction_[target_level+1]->do_restrict(*sp_andra_cu, *sp_fine_nontop);

            //VecScaleAdd(*sp_cu,0.5,*sp_cu,0.5,*sp_andra_cu);
            sp_cu->set_storage_type(sp_fine_nontop->get_storage_type());
            sp_cu->enable_redistribution(sp_fine_nontop->redistribution_enabled());

            return sp_cu;
}

template<typename TDomain, typename TAlgebra>
void SpatialGridTransfer<TDomain, TAlgebra>::set_approx_space(SP_ApproximationSpace approx_space) {
    this->approximation_space_ = approx_space;
}

template<typename TDomain, typename TAlgebra>
void SpatialGridTransfer<TDomain, TAlgebra>::set_prolongation(SP_TransferOperator prolongation) {
    this->prolongation_ = prolongation;
}

template<typename TDomain, typename TAlgebra>
void SpatialGridTransfer<TDomain, TAlgebra>::set_restriction(SP_TransferOperator restriction) {
    this->restriction_ = restriction;
}

template<typename TDomain, typename TAlgebra>
void SpatialGridTransfer<TDomain, TAlgebra>::set_domain(SP_DomainDisc domain_disc) {
    this->domain_disc_ = domain_disc;
}

template<typename TDomain, typename TAlgebra>
void SpatialGridTransfer<TDomain, TAlgebra>::set_transfer(SP_TransferOperator transfer) {
    this->set_prolongation(transfer);
    this->set_restriction(transfer);
}

template<typename TDomain, typename TAlgebra>
void SpatialGridTransfer<TDomain, TAlgebra>::init() {

    const size_t num_level = this->approximation_space_->num_levels();

    if (restriction_ == SPNULL) {
        std::cout << " using StdInjection for restriction " << std::endl;
        //this->restriction_ = make_sp(new StdInjection<TDomain,TAlgebra>(this->approximation_space_));
    }
    {
        //std::cout << " using StdTransfer for restriction " << std::endl;


        SmartPtr<StdTransfer<TDomain,TAlgebra>> stdtransfer  = make_sp(new StdTransfer<TDomain,TAlgebra>());
        //auto dbgWriter = make_sp(new GridFunctionDebugWriter<TDomain,TAlgebra>(this->approximation_space_));
        //stdtransfer->set_debug(dbgWriter);
        //stdtransfer->enable_p1_lagrange_optimization(false);
        stdtransfer->set_restriction_damping(0.25);

        this->restriction_ = stdtransfer;
        //this->andra_restriction_ = stdtransfer;
    }
    if (prolongation_ == SPNULL) {
        std::cout << " using StdTransfer for prolongation " << std::endl;

        //this->prolongation_ = make_sp(new StdInjection<TDomain,TAlgebra>(this->approximation_space_));
        SmartPtr<StdTransfer<TDomain,TAlgebra>> stdtransfer = make_sp(new StdTransfer<TDomain,TAlgebra>());
        //stdtransfer->enable_p1_lagrange_optimization(true);
        //auto dbgWriter = make_sp(new GridFunctionDebugWriter<TDomain,TAlgebra>(this->approximation_space_));
        //stdtransfer->set_debug(dbgWriter);
        this->prolongation_ = stdtransfer;

    }

    const size_t num_const = domain_disc_->num_constraints();
    for( size_t i = 0 ; i < num_const; ++i) {
        this->prolongation_->add_constraint(domain_disc_->constraint(i));
        this->restriction_->add_constraint(domain_disc_->constraint(i));
        //this->andra_restriction_->add_constraint(domain_disc_->constraint(i));
    }

    level_restriction_.resize(num_level);
    //andra_level_restriction_.resize(num_level);
    level_prolongation_.resize(num_level);
    const size_t base_level = 1; // 0 is base level restriction 1 -> 0, prolongation 0 -> 1 both called from lvl 1



    for (size_t level = base_level; level < num_level; ++level) {
        SP_TransferOperator op_restriction = this->restriction_->clone();
        //SP_TransferOperator op_andra_restriction = this->andra_restriction_->clone();
        SP_TransferOperator op_prolongation = this->prolongation_->clone();

        GridLevel coarse_level = GridLevel(level-1, GridLevel::LEVEL, false);
        GridLevel fine_level = GridLevel(level, GridLevel::LEVEL, false);

        op_restriction->set_levels(coarse_level, fine_level);
        //op_andra_restriction->set_levels(coarse_level, fine_level);
        op_prolongation->set_levels(coarse_level, fine_level);

        op_restriction->init();
        //op_andra_restriction->init();
        op_prolongation->init();

        level_restriction_[level] = op_restriction;
        //andra_level_restriction_[level] = op_andra_restriction;
        level_prolongation_[level] = op_prolongation;
    }

    {
    }

}
}
}

#endif