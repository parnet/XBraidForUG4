
template<typename TDomain, typename TAlgebra>
bool SimpleTimeIntegrator<TDomain, TAlgebra>::apply_single_stage
(
	SmartPtr<grid_function_type> u1,
	number t1,
	ConstSmartPtr<grid_function_type> u0,
	number t0
)
{
	LIMEX_PROFILE_FUNC();

	// short-cuts
	GridLevel const &gl = u0->grid_level();
	time_disc_type &tdisc = *tdisc_dep_type::m_spTimeDisc;
	typename base_type::solver_type &solver = *base_type::m_spSolver;

	// create solution vector & right hand side
	SmartPtr<grid_function_type> uold;

	// init solution time series
	SmartPtr<vector_time_series_type> m_spSolTimeSeries;        ///< contains all solutions compute so far
	m_spSolTimeSeries=make_sp(new vector_time_series_type());
	m_spSolTimeSeries->clear();
	m_spSolTimeSeries->push(u0->clone(), t0);

	// init solver (and matrix operator)
	SmartPtr<typename base_type::assembled_operator_type> spAssOp;
	spAssOp = make_sp(new typename base_type::assembled_operator_type(tdisc_dep_type::m_spTimeDisc, gl));
	solver.init(spAssOp);

	// integrate
	double t = t0;
	number currdt = base_type::m_dt;
	int step = 1;

	double final_dt = base_type::m_dt;


	while(!hasTerminated(t, t0, t1)) {

		if (this->debug_writer_valid())
		{
			char debug_name_ext[16]; snprintf(debug_name_ext, 16, "%04d", step);
			this->enter_debug_writer_section(std::string("SimpleTimeIntegrator_step") + debug_name_ext);
		}

		// determine step size
		UG_COND_THROW(currdt < base_type::get_dt_min(), "Time step size below minimum. ABORTING!")
			number dt = std::min(currdt, t1-t);
		final_dt = dt;

		// prepare step
		tdisc.prepare_step(m_spSolTimeSeries, dt);
		if (solver.prepare(*u1) == false)
		{

			currdt *= base_type::get_reduction_factor();
			continue;
		}

		// execute step
		if (solver.apply(*u1))
		{
			//
			// ACCEPT step
			//


			// update time
			t += dt;

			// push updated solution into time series (and continue)
			//SmartPtr<typename base_type::vector_type> utmp = m_spSolTimeSeries->oldest();
			//VecAssign(*utmp, static_cast<typename base_type::vector_type> (*u1) );
			uold = m_spSolTimeSeries->push_discard_oldest(u1->clone(), t).template cast_static<grid_function_type>();
			}
		else
		{
			//
			// REJECT step
			//

			currdt *= base_type::get_reduction_factor();
			continue;
		}

		// consistency check
		if (step == 1 && m_spDerivative.valid())
		{

			*m_spDerivative = *u0;
			m_initial_consistency_error = m_spBanachSpace->distance(*m_spDerivative, *u1);
		}

		if (this->debug_writer_valid())
		{
			this->leave_debug_writer_section();
			char debug_name_ext[16]; snprintf(debug_name_ext, 16, "%04d", step);
			this->write_debug(*u1, std::string("TimeSolution_step") + debug_name_ext);
		}

		step++;
		// tdisc.finish_step_elem(m_spSolTimeSeries, dt);
	}



	if (m_spDerivative.valid())
	{

		VecScaleAdd(static_cast<typename TAlgebra::vector_type&>(*m_spDerivative),
				1.0/final_dt, static_cast<typename TAlgebra::vector_type&>(*u1),
				-1.0/final_dt, static_cast<typename TAlgebra::vector_type&>(*uold));

	}

	m_spSolTimeSeries->clear();

	return true;

};