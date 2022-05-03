/*
 * SWE_Plane_TS_l_rexi_na_sl_nd_etdrk.hpp
 *
 *  Created on: 09 Oct 2017
 *      Author: Pedro Peixoto <pedrosp@ime.usp.br>
 *
 *  Changelog:
 *      based on Martin Schreiber ETD timestepper
 */

#include "../swe_plane_timeintegrators/SWE_Plane_TS_l_rexi_na_sl_nd_etdrk.hpp"

/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Plane_TS_l_rexi_na_sl_nd_etdrk::euler_timestep_update_nonlinear(
		const PlaneData_Spectral &i_h,	///< prognostic variables
		const PlaneData_Spectral &i_u,	///< prognostic variables
		const PlaneData_Spectral &i_v,	///< prognostic variables

		PlaneData_Spectral &o_h_t,	///< time updates
		PlaneData_Spectral &o_u_t,	///< time updates
		PlaneData_Spectral &o_v_t,	///< time updates

		double i_timestamp
)
{
	/*
	 * non-conservative (advective) formulation:
	 *
	 *	h_t = -(u*h)_x - (v*h)_y
	 *	u_t = -g * h_x - u * u_x - v * u_y + f*v
	 *	v_t = -g * h_y - u * v_x - v * v_y - f*u
	 */
	/*
	 * o_h_t = -op.diff_c_x(i_u*i_h) - op.diff_c_y(i_v*i_h);
	 * o_u_t = -i_u*op.diff_c_x(i_u) - i_v*op.diff_c_y(i_u);
	 * o_v_t = -i_u*op.diff_c_x(i_v) - i_v*op.diff_c_y(i_v);
	 */
	// In lagrangian form, the only nonlinearity is the nonlinear divergence
	o_u_t.spectral_set_zero(); //-i_u*op.diff_c_x(i_u) - i_v*op.diff_c_y(i_u);
	o_v_t.spectral_set_zero(); // = 0.0; //-i_u*op.diff_c_x(i_v) - i_v*op.diff_c_y(i_v);

	// linear div only
	if (use_only_linear_divergence)
	{
		o_h_t.spectral_set_zero(); // = 0.0; //-op.diff_c_x(i_u*i_h) - op.diff_c_y(i_v*i_h);
	}
	else
	{
		// nonlinear div
		o_h_t = -i_h*(op.diff_c_x(i_u) + op.diff_c_y(i_v));
		// Smooth spectrum to avoid instability
		if (simVars.misc.use_nonlinear_only_visc != 0)
		{
#if !SWEET_USE_PLANE_SPECTRAL_SPACE
			SWEETError("Implicit diffusion only supported with spectral space activated");
#else
			o_h_t= op.implicit_diffusion(o_h_t, simVars.timecontrol.current_timestep_size*simVars.sim.viscosity, simVars.sim.viscosity_order);
#endif
		}

	}

	if ( simVars.disc.coriolis_treatment == "nonlinear" )
	{
		o_u_t += simVars.sim.plane_rotating_f0 * i_v;
		o_v_t -= simVars.sim.plane_rotating_f0 * i_u;
	}

}



void SWE_Plane_TS_l_rexi_na_sl_nd_etdrk::run_timestep(
		PlaneData_Spectral &io_h,	///< prognostic variables
		PlaneData_Spectral &io_u,	///< prognostic variables
		PlaneData_Spectral &io_v,	///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	if (i_dt <= 0)
		SWEETError("SWE_Plane_TS_l_rexi_na_sl_nd_etdrk: Only constant time step size allowed (Please set --dt)");


	const PlaneDataConfig *planeDataConfig = io_h.planeDataConfig;

	// Tmp vars
	//h, u, v tmp
	PlaneData_Spectral h(planeDataConfig);
	PlaneData_Spectral u(planeDataConfig);
	PlaneData_Spectral v(planeDataConfig);
	//Nonlinear calculation of u,v,h from input
	PlaneData_Spectral FUn_h(planeDataConfig);
	PlaneData_Spectral FUn_u(planeDataConfig);
	PlaneData_Spectral FUn_v(planeDataConfig);


	//Departure points and arrival points
	ScalarDataArray posx_d(io_h.planeDataConfig->physical_array_data_number_of_elements);
	ScalarDataArray posy_d(io_h.planeDataConfig->physical_array_data_number_of_elements);


	Staggering staggering;
	assert(staggering.staggering_type == 'a');

	if (i_simulation_timestamp == 0)
	{
		/*
		 * First time step
		 */
		h_prev = io_h;
		u_prev = io_u;
		v_prev = io_v;
	}


	//Preserve io unmodified
	u = io_u;
	v = io_v;
	h = io_h;

	//////////////////////////////////////////////////////////////
	// Declare quantities common to several timestepping orders //
	//////////////////////////////////////////////////////////////
	// psi1 applied to N(U_0)
	PlaneData_Spectral psi1_FUn_h(planeDataConfig);
	PlaneData_Spectral psi1_FUn_u(planeDataConfig);
	PlaneData_Spectral psi1_FUn_v(planeDataConfig);
	// Interpolated (.5 * psi1NUn + U) to departure points ( )_*
	PlaneData_Spectral Upsi1FUn_h_dep(planeDataConfig);
	PlaneData_Spectral Upsi1FUn_u_dep(planeDataConfig);
	PlaneData_Spectral Upsi1FUn_v_dep(planeDataConfig);





	// Calculate departure points - always force to be second order accurate!
	semiLagrangian.semi_lag_departure_points_settls(
			u_prev.toPhys(),	v_prev.toPhys(),
			u.toPhys(),		v.toPhys(),
			posx_a,		posy_a,
			i_dt,
			posx_d,	posy_d,			// output
			simVars.sim.plane_domain_size,
			&staggering,
			2, //simVars.disc.timestepping_order,

			simVars.disc.semi_lagrangian_max_iterations,
			simVars.disc.semi_lagrangian_convergence_threshold
	);

	////////////////////////////////////
	// Coriolis in advection: f * pos //
	////////////////////////////////////
	// Position x and y
	PlaneData_Spectral fposx_a(planeDataConfig);
	PlaneData_Spectral fposy_a(planeDataConfig);
	PlaneData_Spectral fposx_d(planeDataConfig);
	PlaneData_Spectral fposy_d(planeDataConfig);
	fposx_a.spectral_set_zero();
	fposy_a.spectral_set_zero();
	fposx_d.spectral_set_zero();
	fposy_d.spectral_set_zero();
	if ( simVars.disc.coriolis_treatment == "advection" )
	{
		// The term added (and after subtracted) to (u,v) = f * k x r = (-fy, fx)
		PlaneData_Physical posx_a_phys = Convert_ScalarDataArray_to_PlaneDataPhysical::convert(posx_a, planeDataConfig);
		PlaneData_Physical posy_a_phys = Convert_ScalarDataArray_to_PlaneDataPhysical::convert(posy_a, planeDataConfig);

		PlaneData_Physical fposx_a_phys = simVars.sim.plane_rotating_f0 * posx_a_phys;
		PlaneData_Physical fposy_a_phys = -simVars.sim.plane_rotating_f0 * posy_a_phys;

		fposx_a.loadPlaneDataPhysical(fposx_a_phys);
		fposy_a.loadPlaneDataPhysical(fposy_a_phys);

		// The term added (and after subtracted) to (u,v) = f * k x r = (-fy, fx)
		PlaneData_Physical posx_d_phys = Convert_ScalarDataArray_to_PlaneDataPhysical::convert(posx_d, planeDataConfig);
		PlaneData_Physical posy_d_phys = Convert_ScalarDataArray_to_PlaneDataPhysical::convert(posy_d, planeDataConfig);

		PlaneData_Physical fposx_d_phys = simVars.sim.plane_rotating_f0 * posx_d_phys;
		PlaneData_Physical fposy_d_phys = -simVars.sim.plane_rotating_f0 * posy_d_phys;

		fposx_d.loadPlaneDataPhysical(fposx_d_phys);
		fposy_d.loadPlaneDataPhysical(fposy_d_phys);

	}


	if (timestepping_order == 1 || timestepping_order == 2 || timestepping_order == -2 || timestepping_order == -22)
	{
		/*
		 * U_{1} = \phi_{0}( \Delta t L ) [
		 * 			U_{0}_dep + \Delta t  (\phi_{1}(-\Delta tL) N(U_{0}))_dep.
		 *
		 *\phi_{1}(-\Delta tL)=psi_{1}(\Delta tL)
		 *
		 *F(U)=N(U)
		 *
		 */

		// Calculate term to be interpolated: u+dt*psi_1(dt L)N(U_{0})
		//Calculate N(U_{0})

		if (!use_only_linear_divergence) //Full nonlinear case
		{

			euler_timestep_update_nonlinear(
					h, u, v,
					FUn_h, FUn_u, FUn_v,
					i_simulation_timestamp
			);


			//Apply psi_1 to N(U_{0})
			/////PlaneData_Spectral psi1_FUn_h(planeDataConfig);
			/////PlaneData_Spectral psi1_FUn_u(planeDataConfig);
			/////PlaneData_Spectral psi1_FUn_v(planeDataConfig);

			ts_psi1_rexi.run_timestep(
					FUn_h, FUn_u, FUn_v,
					psi1_FUn_h, psi1_FUn_u, psi1_FUn_v,
					i_dt,
					i_simulation_timestamp
			);


			//Add this to U and interpolate to departure points
			h = h + i_dt*psi1_FUn_h;
			u = u + i_dt*psi1_FUn_u;
			v = v + i_dt*psi1_FUn_v;
		}

		h = sampler2D.bicubic_scalar(h, posx_d, posy_d, -0.5, -0.5);
		u = sampler2D.bicubic_scalar(u, posx_d, posy_d, -0.5, -0.5);
		v = sampler2D.bicubic_scalar(v, posx_d, posy_d, -0.5, -0.5);


		// Add coriolis term to (u,v)
		// IT is known analytically, therefore no need of interpolation // 20197-ifs-documentation-cy47r3-part-iii-dynamics-and-numerical-procedures.pdf
		// Store it in u_with_coriolis, interpolate it and apply phi
		// and let (u,v) as usual to compute nonlinear term in higher orders
		PlaneData_Spectral tmp_u = u;
		if ( simVars.disc.coriolis_treatment == "advection" )
		{
			std::cout << "CCC" << std::endl;
			u = u + fposy_d;
			v = v + fposx_d;

			///u.loadPlaneDataPhysical(u.toPhys() + fposy_d.toPhys());
			///v.loadPlaneDataPhysical(v.toPhys() + fposx_d.toPhys());


			//////////u_with_coriolis = u + fposy;
			//////////v_with_coriolis = v + fposx;

			//////////u_with_coriolis = sampler2D.bicubic_scalar(u_with_coriolis, posx_d, posy_d, -0.5, -0.5);
			//////////v_with_coriolis = sampler2D.bicubic_scalar(v_with_coriolis, posx_d, posy_d, -0.5, -0.5);

			//////////PlaneData_Spectral h_dummy = h;

			////////////Calculate phi_0 of interpolated U
			//////////ts_phi0_rexi.run_timestep(
			//////////		h, u_with_coriolis, v_with_coriolis,
			//////////		h_dummy, u_with_coriolis, v_with_coriolis,
			//////////		i_dt,
			//////////		i_simulation_timestamp
			//////////);

		}


		//Calculate phi_0 of interpolated U
		PlaneData_Spectral phi0_Un_h(planeDataConfig);
		PlaneData_Spectral phi0_Un_u(planeDataConfig);
		PlaneData_Spectral phi0_Un_v(planeDataConfig);
		ts_phi0_rexi.run_timestep(
				h, u, v,
				phi0_Un_h, phi0_Un_u, phi0_Un_v,
				i_dt,
				i_simulation_timestamp
		);

		h = phi0_Un_h;
		u = phi0_Un_u;
		v = phi0_Un_v;




		if ( simVars.disc.coriolis_treatment == "advection" )
		{
			std::cout << "DDD" << std::endl;
			u = u - fposy_d;
			v = v - fposx_d;
			///u.loadPlaneDataPhysical(u.toPhys() - fposy_a.toPhys());
			///v.loadPlaneDataPhysical(v.toPhys() - fposx_a.toPhys());
		}
		std::cout << "AAAA " << (tmp_u - u).toPhys().physical_reduce_max_abs() << std::endl;
	}
	////else
	////{
	////	SWEETError("TODO: This order is not implemented, yet!");
	////}


	// SL-ETD1RK-SETTLS
	if ( timestepping_order == -1 )
	{

		/*
		 * U_{1} = \phi_{0}( \Delta t L ) [
		 * 			U_{0}_dep + \Delta t  (\phi_{1}(-\Delta tL) N(U_{0}))_dep.
		 *
		 *\phi_{1}(-\Delta tL)=psi_{1}(\Delta tL)
		 *
		 *F(U)=N(U)
		 *
		 */

		// Calculate term to be interpolated: u+dt*psi_1(dt L)N(U_{0})
		//Calculate N(U_{0})

		if (!use_only_linear_divergence) //Full nonlinear case
		{

			euler_timestep_update_nonlinear(
					h, u, v,
					FUn_h, FUn_u, FUn_v,
					i_simulation_timestamp
			);


			// Nonlinear term applied to previous time step
			PlaneData_Spectral FUn_h_prev(planeDataConfig);
			PlaneData_Spectral FUn_u_prev(planeDataConfig);
			PlaneData_Spectral FUn_v_prev(planeDataConfig);
			euler_timestep_update_nonlinear(
					h_prev, u_prev, v_prev,
					FUn_h_prev, FUn_u_prev, FUn_v_prev,
					i_simulation_timestamp
			);


			//Apply psi_1 to N(U_{0})
			ts_psi1_rexi.run_timestep(
					2. * FUn_h - FUn_h_prev, 2. * FUn_u - FUn_u_prev, 2. * FUn_v - FUn_v_prev,
					psi1_FUn_h, psi1_FUn_u, psi1_FUn_v,
					i_dt,
					i_simulation_timestamp
			);


			//Add this to U and interpolate to departure points
			h = h + .5 * i_dt*psi1_FUn_h;
			u = u + .5 * i_dt*psi1_FUn_u;
			v = v + .5 * i_dt*psi1_FUn_v;
		}


		h = sampler2D.bicubic_scalar(h, posx_d, posy_d, -0.5, -0.5);
		u = sampler2D.bicubic_scalar(u, posx_d, posy_d, -0.5, -0.5);
		v = sampler2D.bicubic_scalar(v, posx_d, posy_d, -0.5, -0.5);

		// Calculate psi1 of non interpolated nonlinear term
		PlaneData_Spectral psi1_FUn_h_B(planeDataConfig);
		PlaneData_Spectral psi1_FUn_u_B(planeDataConfig);
		PlaneData_Spectral psi1_FUn_v_B(planeDataConfig);
		ts_psi1_rexi.run_timestep(
				FUn_h, FUn_u, FUn_v,
				psi1_FUn_h_B, psi1_FUn_u_B, psi1_FUn_v_B,
				i_dt,
				i_simulation_timestamp
		);

                h = h + .5 * i_dt * psi1_FUn_h_B;
                u = u + .5 * i_dt * psi1_FUn_u_B;
                v = v + .5 * i_dt * psi1_FUn_v_B;

		//Calculate phi_0 of interpolated U
		PlaneData_Spectral phi0_Un_h(planeDataConfig);
		PlaneData_Spectral phi0_Un_u(planeDataConfig);
		PlaneData_Spectral phi0_Un_v(planeDataConfig);
		ts_phi0_rexi.run_timestep(
				h, u, v,
				phi0_Un_h, phi0_Un_u, phi0_Un_v,
				i_dt,
				i_simulation_timestamp
		);



		h = phi0_Un_h;
		u = phi0_Un_u;
		v = phi0_Un_v;

	}


	if ((timestepping_order == -2 || timestepping_order == -22) && !use_only_linear_divergence)
	{

		//Calculate order 1 SL-ETDRK (as above) : {U}_1^{n+1}
		//-----------------------------------------------------
		// Save SL-ETD1RK from above
		PlaneData_Spectral h1(planeDataConfig);
		PlaneData_Spectral u1(planeDataConfig);
		PlaneData_Spectral v1(planeDataConfig);
		h1 = h;
		u1 = u;
		v1 = v;

		/////////////////////////////////Calculate order 1 SL-ETDRK (as above) : {U}_1^{n+1}
		/////////////////////////////////-----------------------------------------------------
		/////////////////////////////////Apply psi_1 to N(U_{0}) (could be stored from the previous if loop)
		///////////////////////////////PlaneData_Spectral psi1_FUn_h(planeDataConfig);
		///////////////////////////////PlaneData_Spectral psi1_FUn_u(planeDataConfig);
		///////////////////////////////PlaneData_Spectral psi1_FUn_v(planeDataConfig);

		///////////////////////////////ts_psi1_rexi.run_timestep(
		///////////////////////////////		FUn_h, FUn_u, FUn_v,
		///////////////////////////////		psi1_FUn_h, psi1_FUn_u, psi1_FUn_v,
		///////////////////////////////		i_dt,
		///////////////////////////////		i_simulation_timestamp
		///////////////////////////////);

		//Add half of this to (original) U and interpolate to departure points
		h = io_h + .5 * i_dt*psi1_FUn_h;
		u = io_u + .5 * i_dt*psi1_FUn_u;
		v = io_v + .5 * i_dt*psi1_FUn_v;

		//Interpolate (.5 * psi1NUn + U) to departure points ( )_*
		//////////////////////////////PlaneData_Spectral Upsi1FUn_h_dep(planeDataConfig);
		//////////////////////////////PlaneData_Spectral Upsi1FUn_u_dep(planeDataConfig);
		//////////////////////////////PlaneData_Spectral Upsi1FUn_v_dep(planeDataConfig);
		Upsi1FUn_h_dep = sampler2D.bicubic_scalar(h, posx_d, posy_d, -0.5, -0.5);
		Upsi1FUn_u_dep = sampler2D.bicubic_scalar(u, posx_d, posy_d, -0.5, -0.5);
		Upsi1FUn_v_dep = sampler2D.bicubic_scalar(v, posx_d, posy_d, -0.5, -0.5);

		//Calculate psi1NU_1
		//-----------------------

		//NU_1
		PlaneData_Spectral FU1_h(planeDataConfig);
		PlaneData_Spectral FU1_u(planeDataConfig);
		PlaneData_Spectral FU1_v(planeDataConfig);
		euler_timestep_update_nonlinear(
				h1, u1, v1,
				FU1_h, FU1_u, FU1_v,
				i_simulation_timestamp
		);

		//Apply psi1
		PlaneData_Spectral psi1_FU1_h(planeDataConfig);
		PlaneData_Spectral psi1_FU1_u(planeDataConfig);
		PlaneData_Spectral psi1_FU1_v(planeDataConfig);
		ts_psi1_rexi.run_timestep(
				FU1_h, FU1_u, FU1_v,
				psi1_FU1_h, psi1_FU1_u, psi1_FU1_v,
				i_dt,
				i_simulation_timestamp
		);

		//Add half of this to already interpolated values
		h = Upsi1FUn_h_dep + .5 * i_dt * psi1_FU1_h;
		u = Upsi1FUn_u_dep + .5 * i_dt * psi1_FU1_u;
		v = Upsi1FUn_v_dep + .5 * i_dt * psi1_FU1_v;

		//Calculate phi_0 of [interpolated (.5 * psi1NUN + U) + .5 * ps1NU1]
		PlaneData_Spectral phi0_Un_h(planeDataConfig);
		PlaneData_Spectral phi0_Un_u(planeDataConfig);
		PlaneData_Spectral phi0_Un_v(planeDataConfig);
		ts_phi0_rexi.run_timestep(
				h, u, v,
				phi0_Un_h, phi0_Un_u, phi0_Un_v,
				i_dt,
				i_simulation_timestamp
		);

		h = phi0_Un_h;
		u = phi0_Un_u;
		v = phi0_Un_v;

	}

	if ( timestepping_order == -22 && !use_only_linear_divergence)
	{

		//Calculate order SL-ETD1RK (as above) : {A}_1^{n+1}
		//-----------------------------------------------------
		// Save SL-ETD1RK from above
		PlaneData_Spectral A_h1(planeDataConfig);
		PlaneData_Spectral A_u1(planeDataConfig);
		PlaneData_Spectral A_v1(planeDataConfig);
		A_h1 = h;
		A_u1 = u;
		A_v1 = v;

		//Calculate order 1 SL-ETDRK (as above) : {U}_1^{n+1}
		//-----------------------------------------------------
		//Apply psi_1 to N(U_{0}) (could be stored from the previous if loop)
		///////////////////////////////PlaneData_Spectral psi1_FUn_h(planeDataConfig);
		///////////////////////////////PlaneData_Spectral psi1_FUn_u(planeDataConfig);
		///////////////////////////////PlaneData_Spectral psi1_FUn_v(planeDataConfig);

		///////////////////////////////ts_psi1_rexi.run_timestep(
		///////////////////////////////		FUn_h, FUn_u, FUn_v,
		///////////////////////////////		psi1_FUn_h, psi1_FUn_u, psi1_FUn_v,
		///////////////////////////////		i_dt,
		///////////////////////////////		i_simulation_timestamp
		///////////////////////////////);

		/////////////////////////////////Add half of this to (original) U and interpolate to departure points
		///////////////////////////////// Also could be stored from previous if loop
		///////////////////////////////h = io_h + .5 * i_dt*psi1_FUn_h;
		///////////////////////////////u = io_u + .5 * i_dt*psi1_FUn_u;
		///////////////////////////////v = io_v + .5 * i_dt*psi1_FUn_v;

		///////////////////////////////// Interpolate (.5 * psi1NUn + U) to departure points ( )_*
		///////////////////////////////// Also could be stored from previous if loop
		///////////////////////////////PlaneData_Spectral Upsi1FUn_h_dep(planeDataConfig);
		///////////////////////////////PlaneData_Spectral Upsi1FUn_u_dep(planeDataConfig);
		///////////////////////////////PlaneData_Spectral Upsi1FUn_v_dep(planeDataConfig);
		///////////////////////////////Upsi1FUn_h_dep = sampler2D.bicubic_scalar(h, posx_d, posy_d, -0.5, -0.5);
		///////////////////////////////Upsi1FUn_u_dep = sampler2D.bicubic_scalar(u, posx_d, posy_d, -0.5, -0.5);
		///////////////////////////////Upsi1FUn_v_dep = sampler2D.bicubic_scalar(v, posx_d, posy_d, -0.5, -0.5);

		//Calculate psi1NA_1
		//-----------------------

		//NU_1
		PlaneData_Spectral FA1_h(planeDataConfig);
		PlaneData_Spectral FA1_u(planeDataConfig);
		PlaneData_Spectral FA1_v(planeDataConfig);
		euler_timestep_update_nonlinear(
				A_h1, A_u1, A_v1,
				FA1_h, FA1_u, FA1_v,
				i_simulation_timestamp
		);

		//Apply psi1
		PlaneData_Spectral psi1_FA1_h(planeDataConfig);
		PlaneData_Spectral psi1_FA1_u(planeDataConfig);
		PlaneData_Spectral psi1_FA1_v(planeDataConfig);
		ts_psi1_rexi.run_timestep(
				FA1_h, FA1_u, FA1_v,
				psi1_FA1_h, psi1_FA1_u, psi1_FA1_v,
				i_dt,
				i_simulation_timestamp
		);

		//Add half of this to already interpolated values
		h = Upsi1FUn_h_dep + .5 * i_dt * psi1_FA1_h;
		u = Upsi1FUn_u_dep + .5 * i_dt * psi1_FA1_u;
		v = Upsi1FUn_v_dep + .5 * i_dt * psi1_FA1_v;

		//Calculate phi_0 of [interpolated (.5 * psi1NUN + U) + .5 * ps1NU1]
		PlaneData_Spectral phi0_Un_h(planeDataConfig);
		PlaneData_Spectral phi0_Un_u(planeDataConfig);
		PlaneData_Spectral phi0_Un_v(planeDataConfig);
		ts_phi0_rexi.run_timestep(
				h, u, v,
				phi0_Un_h, phi0_Un_u, phi0_Un_v,
				i_dt,
				i_simulation_timestamp
		);

		h = phi0_Un_h;
		u = phi0_Un_u;
		v = phi0_Un_v;

	}

	//Aditional steps for SL-ETD2RK (it depends on SL-ETD1RK) - only for full nonlinear case
	if (timestepping_order == 2 && !use_only_linear_divergence)
	{

		/*
		 * U_2^{n+1} =  {U}_1^{n+1}
		 *     + \Delta t\, \varphi_0(\Delta t L)\left[ \psi_2(\Delta t L) N({U}_1^{n+1})
		 *     - \left(\psi_2(\Delta t L)N(U^n) \right)_*^n\right],
		 *
		 *     F(U)=N(U)
		 */

		//Calculate order 1 SL-ETDRK (as above) : {U}_1^{n+1}
		//-----------------------------------------------------
		// Save SL-ETD1RK from above
		PlaneData_Spectral h1(planeDataConfig);
		PlaneData_Spectral u1(planeDataConfig);
		PlaneData_Spectral v1(planeDataConfig);
		h1 = h;
		u1 = u;
		v1 = v;

		//Calculate psi2NU_1
		//-----------------------

		//NU_1
		PlaneData_Spectral FU1_h(planeDataConfig);
		PlaneData_Spectral FU1_u(planeDataConfig);
		PlaneData_Spectral FU1_v(planeDataConfig);
		euler_timestep_update_nonlinear(
				h1, u1, v1,
				FU1_h, FU1_u, FU1_v,
				i_simulation_timestamp
		);


		//Apply psi2
		PlaneData_Spectral psi2_FU1_h(planeDataConfig);
		PlaneData_Spectral psi2_FU1_u(planeDataConfig);
		PlaneData_Spectral psi2_FU1_v(planeDataConfig);
		ts_psi2_rexi.run_timestep(
				FU1_h, FU1_u, FU1_v,
				psi2_FU1_h, psi2_FU1_u, psi2_FU1_v,
				i_dt,
				i_simulation_timestamp
		);

		//Calculate psi2NUn
		PlaneData_Spectral psi2_FUn_h(planeDataConfig);
		PlaneData_Spectral psi2_FUn_u(planeDataConfig);
		PlaneData_Spectral psi2_FUn_v(planeDataConfig);

		ts_psi2_rexi.run_timestep(
				FUn_h, FUn_u, FUn_v,
				psi2_FUn_h, psi2_FUn_u, psi2_FUn_v,
				i_dt,
				i_simulation_timestamp
		);


		//Interpolate psi2NUn to departure points ( )_*

		PlaneData_Spectral psi2FUn_h_dep(planeDataConfig);
		PlaneData_Spectral psi2FUn_u_dep(planeDataConfig);
		PlaneData_Spectral psi2FUn_v_dep(planeDataConfig);
		psi2FUn_h_dep = sampler2D.bicubic_scalar(psi2_FUn_h, posx_d, posy_d, -0.5, -0.5);
		psi2FUn_u_dep = sampler2D.bicubic_scalar(psi2_FUn_u, posx_d, posy_d, -0.5, -0.5);
		psi2FUn_v_dep = sampler2D.bicubic_scalar(psi2_FUn_v, posx_d, posy_d, -0.5, -0.5);


		//psi2NU_1-psi2NUn_dep
		PlaneData_Spectral dif2_h(io_h.planeDataConfig);
		PlaneData_Spectral dif2_u(io_u.planeDataConfig);
		PlaneData_Spectral dif2_v(io_v.planeDataConfig);

		dif2_h = psi2_FU1_h-psi2FUn_h_dep;
		dif2_u = psi2_FU1_u-psi2FUn_u_dep;
		dif2_v = psi2_FU1_v-psi2FUn_v_dep;

		//Apply phi0
		//second order final forcing
		PlaneData_Spectral phi0_dif2_h(planeDataConfig);
		PlaneData_Spectral phi0_dif2_u(planeDataConfig);
		PlaneData_Spectral phi0_dif2_v(planeDataConfig);
		ts_phi0_rexi.run_timestep(
				dif2_h, dif2_u, dif2_v,
				phi0_dif2_h, phi0_dif2_u, phi0_dif2_v,
				i_dt,
				i_simulation_timestamp
		);

		h = h1 + i_dt*phi0_dif2_h;
		u = u1 + i_dt*phi0_dif2_u;
		v = v1 + i_dt*phi0_dif2_v;

	}



	// Save current time step for next step
	h_prev = io_h;
	u_prev = io_u;
	v_prev = io_v;

	io_h = h;
	io_u = u;
	io_v = v;

}



/**
 * Setup
 */
void SWE_Plane_TS_l_rexi_na_sl_nd_etdrk::setup(
	int i_timestepping_order,
	bool i_use_only_linear_divergence
)
{
	timestepping_order = i_timestepping_order;
	use_only_linear_divergence = i_use_only_linear_divergence;

	ts_phi0_rexi.setup(simVars.rexi, "phi0", simVars.timecontrol.current_timestep_size);

	if ((timestepping_order == 1 || timestepping_order == -1 || timestepping_order == -2 || timestepping_order == -22) && !use_only_linear_divergence)
	{
		ts_phi1_rexi.setup(simVars.rexi, "phi1", simVars.timecontrol.current_timestep_size);
		ts_psi1_rexi.setup(simVars.rexi, "psi1", simVars.timecontrol.current_timestep_size);
	}
	else if (timestepping_order == 2 && !use_only_linear_divergence)
	{
		ts_phi1_rexi.setup(simVars.rexi, "phi1", simVars.timecontrol.current_timestep_size);
		ts_phi2_rexi.setup(simVars.rexi, "phi2", simVars.timecontrol.current_timestep_size);

		ts_psi1_rexi.setup(simVars.rexi, "psi1", simVars.timecontrol.current_timestep_size);
		ts_psi2_rexi.setup(simVars.rexi, "psi2", simVars.timecontrol.current_timestep_size);

	}
	else if (timestepping_order == 4 && !use_only_linear_divergence)
	{
		ts_phi1_rexi.setup(simVars.rexi, "phi1", simVars.timecontrol.current_timestep_size*0.5);
		ts_phi2_rexi.setup(simVars.rexi, "phi2", simVars.timecontrol.current_timestep_size*0.5);

		ts_psi1_rexi.setup(simVars.rexi, "psi1", simVars.timecontrol.current_timestep_size*0.5);
		ts_psi2_rexi.setup(simVars.rexi, "psi2", simVars.timecontrol.current_timestep_size*0.5);
		ts_psi3_rexi.setup(simVars.rexi, "psi3", simVars.timecontrol.current_timestep_size*0.5);

		ts_ups0_rexi.setup(simVars.rexi, "phi0", simVars.timecontrol.current_timestep_size);
		ts_ups1_rexi.setup(simVars.rexi, "ups1", simVars.timecontrol.current_timestep_size);
		ts_ups2_rexi.setup(simVars.rexi, "ups2", simVars.timecontrol.current_timestep_size);
		ts_ups3_rexi.setup(simVars.rexi, "ups3", simVars.timecontrol.current_timestep_size);
	}

	if (simVars.disc.space_grid_use_c_staggering)
		SWEETError("SWE_Plane_TS_l_rexi_na_sl_nd_etdrk: Staggering not supported for l_rexi_na_sl_nd_etdrk");

	//with_linear_div_only = i_use_linear_div;

	// Setup sampler for future interpolations
	sampler2D.setup(simVars.sim.plane_domain_size, op.planeDataConfig);

	// Setup semi-lag
	semiLagrangian.setup(simVars.sim.plane_domain_size, op.planeDataConfig);


	PlaneData_Physical tmp_x(op.planeDataConfig);
	tmp_x.physical_update_lambda_array_indices(
			[&](int i, int j, double &io_data)
			{
		io_data = ((double)i)*simVars.sim.plane_domain_size[0]/(double)simVars.disc.space_res_physical[0];
			},
			false
	);

	PlaneData_Physical tmp_y(op.planeDataConfig);
	tmp_y.physical_update_lambda_array_indices(
			[&](int i, int j, double &io_data)
			{
		io_data = ((double)j)*simVars.sim.plane_domain_size[1]/(double)simVars.disc.space_res_physical[1];
			},
			false
	);

	// Initialize arrival points with h position
	ScalarDataArray pos_x = Convert_PlaneDataPhysical_To_ScalarDataArray::physical_convert(tmp_x);
	ScalarDataArray pos_y = Convert_PlaneDataPhysical_To_ScalarDataArray::physical_convert(tmp_y);


	double cell_size_x = simVars.sim.plane_domain_size[0]/(double)simVars.disc.space_res_physical[0];
	double cell_size_y = simVars.sim.plane_domain_size[1]/(double)simVars.disc.space_res_physical[1];

	// Initialize arrival points with h position
	posx_a = pos_x+0.5*cell_size_x;
	posy_a = pos_y+0.5*cell_size_y;

}


SWE_Plane_TS_l_rexi_na_sl_nd_etdrk::SWE_Plane_TS_l_rexi_na_sl_nd_etdrk(
		SimulationVariables &i_simVars,
		PlaneOperators &i_op
)	:
				simVars(i_simVars),
				op(i_op),
				ts_phi0_rexi(simVars, op),
				ts_phi1_rexi(simVars, op),
				ts_phi2_rexi(simVars, op),

				ts_ups0_rexi(simVars, op),
				ts_ups1_rexi(simVars, op),
				ts_ups2_rexi(simVars, op),
				ts_ups3_rexi(simVars, op),

				ts_psi1_rexi(simVars, op),
				ts_psi2_rexi(simVars, op),
				ts_psi3_rexi(simVars, op),

				h_prev(i_op.planeDataConfig),
				u_prev(i_op.planeDataConfig),
				v_prev(i_op.planeDataConfig),

				posx_a(i_op.planeDataConfig->physical_array_data_number_of_elements),
				posy_a(i_op.planeDataConfig->physical_array_data_number_of_elements),

				posx_d(i_op.planeDataConfig->physical_array_data_number_of_elements),
				posy_d(i_op.planeDataConfig->physical_array_data_number_of_elements)
{
}



SWE_Plane_TS_l_rexi_na_sl_nd_etdrk::~SWE_Plane_TS_l_rexi_na_sl_nd_etdrk()
{
}

