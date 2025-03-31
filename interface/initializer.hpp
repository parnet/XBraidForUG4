#ifndef UGPLUGIN_XBRAIDFORUG4_INTERFACE_INITIALIZER_HPP
#define UGPLUGIN_XBRAIDFORUG4_INTERFACE_INITIALIZER_HPP







namespace ug{ namespace xbraid {

    template <typename TDomain, typename TAlgebra>
    class BraidInitializer {
    public:
        //--------------------------------------------------------------------------------------------------------------

        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;

        using T_ApproximationSpace = ApproximationSpace<TDomain> ;
        using SP_ApproximationSpace = SmartPtr<T_ApproximationSpace> ;

        //--------------------------------------------------------------------------------------------------------------

        SP_GridFunction m_u0 = SPNULL; // u0 and t_start for star value problem
        double t_start = 0.0;

        //--------------------------------------------------------------------------------------------------------------
    protected:
        BraidInitializer() = default;
    public:
        virtual ~BraidInitializer() = default;

        //--------------------------------------------------------------------------------------------------------------
        virtual void initialize(SP_GridFunction& u, number time) = 0;

        virtual void set_start_values(SP_GridFunction u, number time) {


#ifdef FEATURE_SPATIAL_REFINE
        SP_ApproximationSpace approx_space = u->approx_space();

        size_t num_level = approx_space->num_levels();
        std::cout << "num-level = " << num_level << std::endl;
        std::cout << "target --> " << num_level -1<< std::endl;
        SP_GridFunction grid_function = make_sp(new T_GridFunction(approx_space,num_level -1,false));
        std::cout << "gf size =   " << grid_function->size();
        std::cout << "u  size =   " << u->size();

        const size_t usize = u->size();
        for( size_t i = 0; i < usize; ++i) {
            (*grid_function)[i] = (*u)[i];
        }
        grid_function->enable_redistribution(u->redistribution_enabled());
        grid_function->set_storage_type(u->get_storage_mask());
        grid_function->set_layouts(u->layouts());

        this->m_u0 = grid_function;
#else
            this->m_u0 = u;
#endif

            this->t_start = time;
        }


        //--------------------------------------------------------------------------------------------------------------
    };

}}

#endif