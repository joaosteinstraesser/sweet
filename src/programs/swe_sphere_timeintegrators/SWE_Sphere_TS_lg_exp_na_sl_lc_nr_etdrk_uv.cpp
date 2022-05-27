/*
 * SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etdrk_uv
 *
 * Created on: 24 Mar 2022
 *     Author: Joao Steinstraesser <joao.steinstraesser@usp.br>
 *
 */

#include "SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etdrk_uv.hpp"


bool SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etdrk_uv::implements_timestepping_method(const std::string &i_timestepping_method)
{
#if 0
	/*
	 * Should contain _exp and _settls as well as _uv to indicate vorticity-divergence formulation
	 */
	return (
		(i_timestepping_method.find("_exp") != std::string::npos)		&&
		(i_timestepping_method.find("_etd") != std::string::npos)	&&
		(i_timestepping_method.find("_uv") != std::string::npos)		&&
		true
	);
#endif

	if (	i_timestepping_method == "lg_exp_na_sl_lc_nr_etdrk_uv"	||
			i_timestepping_method == "lg_exp_na_sl_lc_etdrk_uv"	||
			false
	)
		return true;

	return false;
}


std::string SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etdrk_uv::string_id()
{
	return "lg_exp_na_sl_lc_nr_etdrk_uv";
}


void SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etdrk_uv::setup_auto()
{
	if (simVars.sim.sphere_use_fsphere)
		SWEETError("TODO: Not yet supported");

	NLRemainderTreatment_enum nonlinear_remainder_treatment = NLRemainderTreatment_enum::NL_REMAINDER_NONLINEAR;

	if (simVars.disc.timestepping_method == "lg_exp_na_sl_lc_nr_etdrk_uv")
	{
		nonlinear_remainder_treatment = NLRemainderTreatment_enum::NL_REMAINDER_NONLINEAR;
	}
	else if (simVars.disc.timestepping_method == "lg_exp_na_sl_lc_etd_uv")
	{
		nonlinear_remainder_treatment = NLRemainderTreatment_enum::NL_REMAINDER_IGNORE;
	}
	else
	{
		SWEETError(std::string("Timestepping method '")+simVars.disc.timestepping_method+std::string("' not known"));
	}

	setup(
			simVars.rexi,
			simVars.disc.timestepping_order,
			simVars.disc.timestepping_order2,
			simVars.timecontrol.current_timestep_size,

			nonlinear_remainder_treatment
		);
}



void SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etdrk_uv::print_help()
{
	std::cout << "	Exponential SL ETD:" << std::endl;
	//std::cout << "		+ lg_exp_na_sl_etd_uv" << std::endl;
	//std::cout << "		+ lg_exp_na_sl_nr_etdrk_uv" << std::endl;
	//std::cout << "		+ lg_exp_na_sl_lc_etd_uv" << std::endl;
	std::cout << "		+ lg_exp_na_sl_lc_nr_etdrk_uv" << std::endl;
#if 0
	std::cout << "		+ l_exp_na_sl_etd_uv" << std::endl;
	std::cout << "		+ l_exp_na_sl_nr_etdrk_uv" << std::endl;
#endif
}



void SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etdrk_uv::run_timestep(
		SphereData_Spectral &io_U_phi,	///< prognostic variables
		SphereData_Spectral &io_U_vrt,	///< prognostic variables
		SphereData_Spectral &io_U_div,	///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	const SphereData_Config *sphereDataConfig = io_U_phi.sphereDataConfig;

////	if (timestepping_order != 2)
////	{
////		SWEETError("Only 2nd order supported");
////	}

	///const SphereData_Spectral &U_phi = io_U_phi;
	///const SphereData_Spectral &U_vrt = io_U_vrt;
	///const SphereData_Spectral &U_div = io_U_div;

	// Keep io data unchanged
	SphereData_Spectral U_phi = io_U_phi;
	SphereData_Spectral U_vrt = io_U_vrt;
	SphereData_Spectral U_div = io_U_div;

	if (i_simulation_timestamp == 0)
	{
		/*
		 * First time step:
		 * Simply backup existing fields for multi-step parts of this algorithm.
		 */
		U_phi_prev = U_phi;
		U_vrt_prev = U_vrt;
		U_div_prev = U_div;
	}


	/*************************************************************************************************
	 * Step 1) Compute departure points
	 *************************************************************************************************/
	SphereData_Physical U_u_lon_prev, U_v_lat_prev;
	ops.vrtdiv_to_uv(U_vrt_prev, U_div_prev, U_u_lon_prev, U_v_lat_prev);

	SphereData_Physical U_u_lon, U_v_lat;
	ops.vrtdiv_to_uv(U_vrt, U_div, U_u_lon, U_v_lat);

	double dt_div_radius = i_dt / simVars.sim.sphere_radius;

	// Calculate departure points
	ScalarDataArray pos_lon_d, pos_lat_d;
	semiLagrangian.semi_lag_departure_points_settls_specialized(
			dt_div_radius*U_u_lon_prev, dt_div_radius*U_v_lat_prev,
			dt_div_radius*U_u_lon, dt_div_radius*U_v_lat,
			pos_lon_d, pos_lat_d		// OUTPUT
		);


	// N(Un) = Nonlinear term applied to Un
	SphereData_Spectral FUn_phi(sphereDataConfig, 0);
	SphereData_Spectral FUn_vrt(sphereDataConfig, 0);
	SphereData_Spectral FUn_div(sphereDataConfig, 0);

	// psi1 applied to N(Un)
	SphereData_Spectral psi1_FUn_phi(sphereDataConfig, 0);
	SphereData_Spectral psi1_FUn_vrt(sphereDataConfig, 0);
	SphereData_Spectral psi1_FUn_div(sphereDataConfig, 0);

	// Interpolated (.5 * psi1NUn + U) to departure points ( )_* //
	SphereData_Spectral psi1_FUn_phi_D(sphereDataConfig, 0);
	SphereData_Spectral psi1_FUn_vrt_D(sphereDataConfig, 0);
	SphereData_Spectral psi1_FUn_div_D(sphereDataConfig, 0);



	// All the schemes use the result provided by ETD1RK
	if (timestepping_order == 1 || timestepping_order == 2 || timestepping_order == -2 || timestepping_order == -22 || timestepping_order == 4)
	{

int test = 1;

///ops.fg = ops.fg * 0;
		////////////////////
		// Compute N(U_n) //
		////////////////////
		ts_ln_erk_split_uv.euler_timestep_update_lc(
				U_phi, U_vrt, U_div,
				FUn_phi, FUn_vrt, FUn_div,
				i_simulation_timestamp
		);

///////////std::cout << "PHI" << std::endl;
///////////FUn_phi.toPhys().physical_print();
///////////std::cout << "VRT" << std::endl;
///////////FUn_vrt.toPhys().physical_print();
///////////std::cout << "DIV" << std::endl;
///////////FUn_div.toPhys().physical_print();
///////////std::cout << "END" << std::endl;

		if (nonlinear_remainder_treatment == NLRemainderTreatment_enum::NL_REMAINDER_NONLINEAR)
		{
			ts_ln_erk_split_uv.euler_timestep_update_nr(
					U_phi, U_vrt, U_div,
					FUn_phi, FUn_vrt, FUn_div,
					i_simulation_timestamp
			);
		}



///ops.fg.physical_print();


if (test == 1)
{


		///////////////////////////
		// Apply psi_1 to N(U_n) //
		///////////////////////////
		ts_psi1_exp.run_timestep(
				FUn_phi, FUn_vrt, FUn_div,
				psi1_FUn_phi, psi1_FUn_vrt, psi1_FUn_div,
				i_dt,
				i_simulation_timestamp
			);


		////////////////////////////////////////////////////////
		// Add this to Un and interpolate to departure points //
		////////////////////////////////////////////////////////
		U_phi = U_phi + i_dt * psi1_FUn_phi;
		U_vrt = U_vrt + i_dt * psi1_FUn_vrt;
		U_div = U_div + i_dt * psi1_FUn_div;

		SphereData_Spectral U_phi_D, U_vrt_D, U_div_D;
		semiLagrangian.apply_sl_timeintegration_uv(
				ops,
				U_phi, U_vrt, U_div,
				pos_lon_d, pos_lat_d,
				U_phi_D, U_vrt_D, U_div_D
			);


		//////////////////////////////////////////
		// Apply phi_0 to interpolated solution //
		//////////////////////////////////////////
		SphereData_Spectral phi0_Un_phi(sphereDataConfig, 0);
		SphereData_Spectral phi0_Un_vrt(sphereDataConfig, 0);
		SphereData_Spectral phi0_Un_div(sphereDataConfig, 0);

		ts_phi0_exp.run_timestep(
				U_phi_D, U_vrt_D, U_div_D,
				phi0_Un_phi, phi0_Un_vrt, phi0_Un_div,
				i_dt,
				i_simulation_timestamp
			);

		U_phi = phi0_Un_phi;
		U_vrt = phi0_Un_vrt;
		U_div = phi0_Un_div;
}
else if (test == 2)
{

		///////////////////////////
		// Apply psi_1 to N(U_n) //
		///////////////////////////
		ts_psi1_exp.run_timestep(
				FUn_phi, FUn_vrt, FUn_div,
				psi1_FUn_phi, psi1_FUn_vrt, psi1_FUn_div,
				i_dt,
				i_simulation_timestamp
			);

		//////////////////////////////////////////
		// Apply phi_0 to nonlinear term        //
		//////////////////////////////////////////
		SphereData_Spectral phi0_FUn_phi(sphereDataConfig, 0);
		SphereData_Spectral phi0_FUn_vrt(sphereDataConfig, 0);
		SphereData_Spectral phi0_FUn_div(sphereDataConfig, 0);
		ts_phi0_exp.run_timestep(
				i_dt * psi1_FUn_phi, i_dt * psi1_FUn_vrt, i_dt * psi1_FUn_div,
				phi0_FUn_phi, phi0_FUn_vrt, phi0_FUn_div,
				i_dt,
				i_simulation_timestamp
			);

		// Interpolate nonlinear term
		SphereData_Spectral phi0_FUn_phi_D(sphereDataConfig, 0);
		SphereData_Spectral phi0_FUn_vrt_D(sphereDataConfig, 0);
		SphereData_Spectral phi0_FUn_div_D(sphereDataConfig, 0);
		semiLagrangian.apply_sl_timeintegration_uv(
				ops,
				phi0_FUn_phi, phi0_FUn_vrt, phi0_FUn_div,
				pos_lon_d, pos_lat_d,
				phi0_FUn_phi_D, phi0_FUn_vrt_D, phi0_FUn_div_D
			);


		// Interpolated solution
		SphereData_Spectral U_phi_D, U_vrt_D, U_div_D;
		semiLagrangian.apply_sl_timeintegration_uv(
				ops,
				U_phi, U_vrt, U_div,
				pos_lon_d, pos_lat_d,
				U_phi_D, U_vrt_D, U_div_D
			);

		//////////////////////////////////////////
		// Apply phi_0 to interpolated solution //
		//////////////////////////////////////////
		SphereData_Spectral phi0_Un_phi(sphereDataConfig, 0);
		SphereData_Spectral phi0_Un_vrt(sphereDataConfig, 0);
		SphereData_Spectral phi0_Un_div(sphereDataConfig, 0);

		ts_phi0_exp.run_timestep(
				U_phi_D, U_vrt_D, U_div_D,
				phi0_Un_phi, phi0_Un_vrt, phi0_Un_div,
				i_dt,
				i_simulation_timestamp
			);

		////////////////////////////////////////////////////////
		// Sum up everything //
		////////////////////////////////////////////////////////
		U_phi = phi0_Un_phi + phi0_FUn_phi_D;
		U_vrt = phi0_Un_vrt + phi0_FUn_vrt_D;
		U_div = phi0_Un_div + phi0_FUn_div_D;


}
else if (test == 3)
{

		///////////////////////////
		// Apply psi_1 to N(U_n) //
		///////////////////////////
		ts_psi1_exp.run_timestep(
				FUn_phi, FUn_vrt, FUn_div,
				psi1_FUn_phi, psi1_FUn_vrt, psi1_FUn_div,
				i_dt,
				i_simulation_timestamp
			);

		// Interpolate nonlinear term
		SphereData_Spectral psi1_FUn_phi_D(sphereDataConfig, 0);
		SphereData_Spectral psi1_FUn_vrt_D(sphereDataConfig, 0);
		SphereData_Spectral psi1_FUn_div_D(sphereDataConfig, 0);
		semiLagrangian.apply_sl_timeintegration_uv(
				ops,
				psi1_FUn_phi, psi1_FUn_vrt, psi1_FUn_div,
				pos_lon_d, pos_lat_d,
				psi1_FUn_phi_D, psi1_FUn_vrt_D, psi1_FUn_div_D
			);

		// Aplly psi0 to interpolated nonlinear term
		SphereData_Spectral phi0_FUn_phi(sphereDataConfig, 0);
		SphereData_Spectral phi0_FUn_vrt(sphereDataConfig, 0);
		SphereData_Spectral phi0_FUn_div(sphereDataConfig, 0);
		ts_phi0_exp.run_timestep(
				i_dt * psi1_FUn_phi_D, i_dt * psi1_FUn_vrt_D, i_dt * psi1_FUn_div_D,
				phi0_FUn_phi, phi0_FUn_vrt, phi0_FUn_div,
				i_dt,
				i_simulation_timestamp
			);


		// Apply phi0 to solution
		SphereData_Spectral phi0_Un_phi(sphereDataConfig, 0);
		SphereData_Spectral phi0_Un_vrt(sphereDataConfig, 0);
		SphereData_Spectral phi0_Un_div(sphereDataConfig, 0);
		ts_phi0_exp.run_timestep(
				U_phi, U_vrt, U_div,
				phi0_Un_phi, phi0_Un_vrt, phi0_Un_div,
				i_dt,
				i_simulation_timestamp
			);


		// Interpolate solution
		SphereData_Spectral phi0_Un_phi_D(sphereDataConfig, 0);
		SphereData_Spectral phi0_Un_vrt_D(sphereDataConfig, 0);
		SphereData_Spectral phi0_Un_div_D(sphereDataConfig, 0);
		semiLagrangian.apply_sl_timeintegration_uv(
				ops,
				phi0_Un_phi, phi0_Un_vrt, phi0_Un_div,
				pos_lon_d, pos_lat_d,
				phi0_Un_phi_D, phi0_Un_vrt_D, phi0_Un_div_D
			);


		//////////////////////////////////////////
		// Sum up everything //
		//////////////////////////////////////////
		U_phi = phi0_Un_phi_D + phi0_FUn_phi;
		U_vrt = phi0_Un_vrt_D + phi0_FUn_vrt;
		U_div = phi0_Un_div_D + phi0_FUn_div;


}
else if (test == 4)
{

		///////////////////////////
		// Apply psi_1 to N(U_n) //
		///////////////////////////
		ts_psi1_exp.run_timestep(
				FUn_phi, FUn_vrt, FUn_div,
				psi1_FUn_phi, psi1_FUn_vrt, psi1_FUn_div,
				i_dt,
				i_simulation_timestamp
			);


		////////////////////////////////////////////////////////
		// Add this to Un                                     //
		////////////////////////////////////////////////////////
		U_phi = U_phi + i_dt * psi1_FUn_phi;
		U_vrt = U_vrt + i_dt * psi1_FUn_vrt;
		U_div = U_div + i_dt * psi1_FUn_div;

		// Apply phi0 to everything
		SphereData_Spectral phi0_Un_phi(sphereDataConfig, 0);
		SphereData_Spectral phi0_Un_vrt(sphereDataConfig, 0);
		SphereData_Spectral phi0_Un_div(sphereDataConfig, 0);
		ts_phi0_exp.run_timestep(
				U_phi, U_vrt, U_div,
				phi0_Un_phi, phi0_Un_vrt, phi0_Un_div,
				i_dt,
				i_simulation_timestamp
			);

		// Interpolate
		SphereData_Spectral U_phi_D, U_vrt_D, U_div_D;
		semiLagrangian.apply_sl_timeintegration_uv(
				ops,
				U_phi, U_vrt, U_div,
				pos_lon_d, pos_lat_d,
				U_phi_D, U_vrt_D, U_div_D
			);

		U_phi = U_phi_D;
		U_vrt = U_vrt_D;
		U_div = U_div_D;


}
else if (test == 5)
{

		// Interpolate nonlinear term
		SphereData_Spectral FUn_phi_D(sphereDataConfig, 0);
		SphereData_Spectral FUn_vrt_D(sphereDataConfig, 0);
		SphereData_Spectral FUn_div_D(sphereDataConfig, 0);
		semiLagrangian.apply_sl_timeintegration_uv(
				ops,
				FUn_phi, FUn_vrt, FUn_div,
				pos_lon_d, pos_lat_d,
				FUn_phi_D, FUn_vrt_D, FUn_div_D
			);

		///////////////////////////
		// Apply psi_1 to N(U_n) //
		///////////////////////////
		ts_psi1_exp.run_timestep(
				FUn_phi_D, FUn_vrt_D, FUn_div_D,
				psi1_FUn_phi, psi1_FUn_vrt, psi1_FUn_div,
				i_dt,
				i_simulation_timestamp
			);


		// Interpolate solution
		SphereData_Spectral U_phi_D, U_vrt_D, U_div_D;
		semiLagrangian.apply_sl_timeintegration_uv(
				ops,
				U_phi, U_vrt, U_div,
				pos_lon_d, pos_lat_d,
				U_phi_D, U_vrt_D, U_div_D
			);


		////////////////////////////////////////////////////////
		// Sum up everything //
		////////////////////////////////////////////////////////
		U_phi = U_phi_D + i_dt * psi1_FUn_phi_D;
		U_vrt = U_vrt_D + i_dt * psi1_FUn_vrt_D;
		U_div = U_div_D + i_dt * psi1_FUn_div_D;


		//////////////////////////////////////////
		// Apply phi_0 //
		//////////////////////////////////////////
		SphereData_Spectral phi0_Un_phi(sphereDataConfig, 0);
		SphereData_Spectral phi0_Un_vrt(sphereDataConfig, 0);
		SphereData_Spectral phi0_Un_div(sphereDataConfig, 0);

		ts_phi0_exp.run_timestep(
				U_phi, U_vrt, U_div,
				phi0_Un_phi, phi0_Un_vrt, phi0_Un_div,
				i_dt,
				i_simulation_timestamp
			);

		U_phi = phi0_Un_phi;
		U_vrt = phi0_Un_vrt;
		U_div = phi0_Un_div;


}

else if (test == 6)
{

int p = this->nb_subintegrals;


		////////////////////////////////////////////////////////
		// Interpolate Un //
		////////////////////////////////////////////////////////
		U_phi = U_phi /*+ i_dt * psi1_FUn_phi*/;
		U_vrt = U_vrt /*+ i_dt * psi1_FUn_vrt*/;
		U_div = U_div /*+ i_dt * psi1_FUn_div*/;

		//SphereData_Spectral U_phi_D, U_vrt_D, U_div_D;
		SphereData_Spectral U_phi_D(sphereDataConfig);
		SphereData_Spectral U_vrt_D(sphereDataConfig);
		SphereData_Spectral U_div_D(sphereDataConfig);
		U_phi_D.spectral_set_zero();
		U_vrt_D.spectral_set_zero();
		U_div_D.spectral_set_zero();
		semiLagrangian.apply_sl_timeintegration_uv(
				ops,
				U_phi, U_vrt, U_div,
				pos_lon_d, pos_lat_d,
				U_phi_D, U_vrt_D, U_div_D
			);

		////////////////////////////////////////////////////////
		// Interpolate FUn //
		////////////////////////////////////////////////////////
		//SphereData_Spectral FUn_phi_D, FUn_vrt_D, FUn_div_D;
		SphereData_Spectral FUn_phi_D(sphereDataConfig);
		SphereData_Spectral FUn_vrt_D(sphereDataConfig);
		SphereData_Spectral FUn_div_D(sphereDataConfig);
		FUn_phi_D.spectral_set_zero();
		FUn_vrt_D.spectral_set_zero();
		FUn_div_D.spectral_set_zero();

		for (int i = 0; i < p; i++)
		{

			double r = (double)i / (double)p;
			double r1 = (double)(i + 1.) / (double)p;
			double dti = r * i_dt;
			double dti1 = r1 * i_dt;

			SphereData_Spectral psi1_FUi_phi(sphereDataConfig);
			SphereData_Spectral psi1_FUi_vrt(sphereDataConfig);
			SphereData_Spectral psi1_FUi_div(sphereDataConfig);
			SphereData_Spectral psi1_FUi1_phi(sphereDataConfig);
			SphereData_Spectral psi1_FUi1_vrt(sphereDataConfig);
			SphereData_Spectral psi1_FUi1_div(sphereDataConfig);
			///////////////////////////
			// Apply psi_1 to N(U_n) //
			///////////////////////////
			ts_psi1_exp.run_timestep(
					FUn_phi, FUn_vrt, FUn_div,
					psi1_FUi_phi, psi1_FUi_vrt, psi1_FUi_div,
					dti,
					i_simulation_timestamp
				);
			ts_psi1_exp.run_timestep(
					FUn_phi, FUn_vrt, FUn_div,
					psi1_FUi1_phi, psi1_FUi1_vrt, psi1_FUi1_div,
					dti1,
					i_simulation_timestamp
				);
			psi1_FUn_phi = dti1 * psi1_FUi1_phi - dti * psi1_FUi_phi;
			psi1_FUn_vrt = dti1 * psi1_FUi1_vrt - dti * psi1_FUi_vrt;
			psi1_FUn_div = dti1 * psi1_FUi1_div - dti * psi1_FUi_div;


			//////////SphereData_Physical U_u_lon_i, U_v_lat_i;
			///////////////ops.vrtdiv_to_uv(U_vrt, U_div, U_u_lon, U_v_lat);
			//////////U_u_lon_i = (1. + r) * U_u_lon - r * U_u_lon_prev;
			//////////U_v_lat_i = (1. + r) * U_v_lat - r * U_v_lat_prev;


			//////////double dt_div_radius_aux = (1. - r) * i_dt / simVars.sim.sphere_radius;

			//////////std::cout << "Computing departure points for i = " << i << "/" << p << std::endl;
			//////////// Calculate departure points
			//////////ScalarDataArray pos_lon_d_2, pos_lat_d_2;
			//////////semiLagrangian.semi_lag_departure_points_settls_specialized(
			//////////		dt_div_radius_aux*U_u_lon_prev, dt_div_radius_aux*U_v_lat_prev,
			//////////		dt_div_radius_aux*U_u_lon_i, dt_div_radius_aux*U_v_lat_i,
			//////////		pos_lon_d_2, pos_lat_d_2		// OUTPUT
			//////////	);

			ScalarDataArray pos_lon_d_aux, pos_lat_d_aux;
			pos_lon_d_aux = r * semiLagrangian.pos_lon_A + (1. - r) * pos_lon_d;
			pos_lat_d_aux = r * semiLagrangian.pos_lat_A + (1. - r) * pos_lat_d;


			std::cout << "Applying SL for p = " << p << std::endl;
			SphereData_Spectral FUn_phi_D_aux, FUn_vrt_D_aux, FUn_div_D_aux;
			semiLagrangian.apply_sl_timeintegration_uv(
					ops,
					psi1_FUn_phi, psi1_FUn_vrt, psi1_FUn_div,
					pos_lon_d_aux, pos_lat_d_aux,
					FUn_phi_D_aux, FUn_vrt_D_aux, FUn_div_D_aux
				);

			std::cout << "Updating nonlinear term for p = " << p << std::endl;
			FUn_phi_D = FUn_phi_D +  FUn_phi_D_aux;
			FUn_vrt_D = FUn_vrt_D +  FUn_vrt_D_aux;
			FUn_div_D = FUn_div_D +  FUn_div_D_aux;

		}

		// Add this to to interpolated solution
		U_phi_D = U_phi_D + FUn_phi_D;
		U_vrt_D = U_vrt_D + FUn_vrt_D;
		U_div_D = U_div_D + FUn_div_D;


		//////////////////////////////////////////
		// Apply phi_0 to interpolated solution //
		//////////////////////////////////////////
		SphereData_Spectral phi0_Un_phi(sphereDataConfig, 0);
		SphereData_Spectral phi0_Un_vrt(sphereDataConfig, 0);
		SphereData_Spectral phi0_Un_div(sphereDataConfig, 0);

		ts_phi0_exp.run_timestep(
				U_phi_D, U_vrt_D, U_div_D,
				phi0_Un_phi, phi0_Un_vrt, phi0_Un_div,
				i_dt,
				i_simulation_timestamp
			);

		U_phi = phi0_Un_phi;
		U_vrt = phi0_Un_vrt;
		U_div = phi0_Un_div;




}

else if (test == 7)
{

int p = this->nb_subintegrals;


		////////////////////////////////////////////////////////
		// Interpolate Un //
		////////////////////////////////////////////////////////
		U_phi = U_phi /*+ i_dt * psi1_FUn_phi*/;
		U_vrt = U_vrt /*+ i_dt * psi1_FUn_vrt*/;
		U_div = U_div /*+ i_dt * psi1_FUn_div*/;

		//SphereData_Spectral U_phi_D, U_vrt_D, U_div_D;
		SphereData_Spectral U_phi_D(sphereDataConfig);
		SphereData_Spectral U_vrt_D(sphereDataConfig);
		SphereData_Spectral U_div_D(sphereDataConfig);
		U_phi_D.spectral_set_zero();
		U_vrt_D.spectral_set_zero();
		U_div_D.spectral_set_zero();
		semiLagrangian.apply_sl_timeintegration_uv(
				ops,
				U_phi, U_vrt, U_div,
				pos_lon_d, pos_lat_d,
				U_phi_D, U_vrt_D, U_div_D
			);

		////////////////////////////////////////////////////////
		// Interpolate FUn //
		////////////////////////////////////////////////////////
		//SphereData_Spectral FUn_phi_D, FUn_vrt_D, FUn_div_D;
		SphereData_Spectral FUn_phi_D(sphereDataConfig);
		SphereData_Spectral FUn_vrt_D(sphereDataConfig);
		SphereData_Spectral FUn_div_D(sphereDataConfig);
		FUn_phi_D.spectral_set_zero();
		FUn_vrt_D.spectral_set_zero();
		FUn_div_D.spectral_set_zero();

		for (int i = 0; i < p; i++)
		{

			double r = (double)i / (double)p;
			double r1 = (double)(i + 1.) / (double)p;
			double dti = r * i_dt;
			double dti1 = r1 * i_dt;

			///SphereData_Spectral psi1_FUi_phi(sphereDataConfig);
			///SphereData_Spectral psi1_FUi_vrt(sphereDataConfig);
			///SphereData_Spectral psi1_FUi_div(sphereDataConfig);
			SphereData_Spectral psi1_FUi1_phi(sphereDataConfig);
			SphereData_Spectral psi1_FUi1_vrt(sphereDataConfig);
			SphereData_Spectral psi1_FUi1_div(sphereDataConfig);
			///////////////////////////
			// Apply psi_1 to N(U_n) //
			///////////////////////////
			///ts_psi1_exp.run_timestep(
			///		FUn_phi, FUn_vrt, FUn_div,
			///		psi1_FUi_phi, psi1_FUi_vrt, psi1_FUi_div,
			///		dti,
			///		i_simulation_timestamp
			///	);
			ts_psi1_exp.run_timestep(
					FUn_phi, FUn_vrt, FUn_div,
					psi1_FUi1_phi, psi1_FUi1_vrt, psi1_FUi1_div,
					dti1,
					i_simulation_timestamp
				);
			psi1_FUn_phi =  1. / (double)p * i_dt * psi1_FUi1_phi;
			psi1_FUn_vrt =  1. / (double)p * i_dt * psi1_FUi1_vrt;
			psi1_FUn_div =  1. / (double)p * i_dt * psi1_FUi1_div;


			//////////SphereData_Physical U_u_lon_i, U_v_lat_i;
			///////////////ops.vrtdiv_to_uv(U_vrt, U_div, U_u_lon, U_v_lat);
			//////////U_u_lon_i = (1. + r) * U_u_lon - r * U_u_lon_prev;
			//////////U_v_lat_i = (1. + r) * U_v_lat - r * U_v_lat_prev;


			//////////double dt_div_radius_aux = (1. - r) * i_dt / simVars.sim.sphere_radius;

			//////////std::cout << "Computing departure points for i = " << i << "/" << p << std::endl;
			//////////// Calculate departure points
			//////////ScalarDataArray pos_lon_d_2, pos_lat_d_2;
			//////////semiLagrangian.semi_lag_departure_points_settls_specialized(
			//////////		dt_div_radius_aux*U_u_lon_prev, dt_div_radius_aux*U_v_lat_prev,
			//////////		dt_div_radius_aux*U_u_lon_i, dt_div_radius_aux*U_v_lat_i,
			//////////		pos_lon_d_2, pos_lat_d_2		// OUTPUT
			//////////	);

			ScalarDataArray pos_lon_d_aux, pos_lat_d_aux;
			pos_lon_d_aux = r * semiLagrangian.pos_lon_A + (1. - r) * pos_lon_d;
			pos_lat_d_aux = r * semiLagrangian.pos_lat_A + (1. - r) * pos_lat_d;
			///pos_lon_d_aux =  pos_lon_d;
			///pos_lat_d_aux =  pos_lat_d;


			std::cout << "Applying SL for p = " << p << std::endl;
			SphereData_Spectral FUn_phi_D_aux, FUn_vrt_D_aux, FUn_div_D_aux;
			semiLagrangian.apply_sl_timeintegration_uv(
					ops,
					psi1_FUn_phi, psi1_FUn_vrt, psi1_FUn_div,
					pos_lon_d_aux, pos_lat_d_aux,
					FUn_phi_D_aux, FUn_vrt_D_aux, FUn_div_D_aux
				);

			std::cout << "Updating nonlinear term for p = " << p << std::endl;
			FUn_phi_D = FUn_phi_D +  FUn_phi_D_aux;
			FUn_vrt_D = FUn_vrt_D +  FUn_vrt_D_aux;
			FUn_div_D = FUn_div_D +  FUn_div_D_aux;

		}

		// Add this to to interpolated solution
		U_phi_D = U_phi_D + FUn_phi_D;
		U_vrt_D = U_vrt_D + FUn_vrt_D;
		U_div_D = U_div_D + FUn_div_D;


		//////////////////////////////////////////
		// Apply phi_0 to interpolated solution //
		//////////////////////////////////////////
		SphereData_Spectral phi0_Un_phi(sphereDataConfig, 0);
		SphereData_Spectral phi0_Un_vrt(sphereDataConfig, 0);
		SphereData_Spectral phi0_Un_div(sphereDataConfig, 0);

		ts_phi0_exp.run_timestep(
				U_phi_D, U_vrt_D, U_div_D,
				phi0_Un_phi, phi0_Un_vrt, phi0_Un_div,
				i_dt,
				i_simulation_timestamp
			);

		U_phi = phi0_Un_phi;
		U_vrt = phi0_Un_vrt;
		U_div = phi0_Un_div;




}




		U_phi.toPhys().physical_print();

	}
	// SL-ETD1RK-bis (2nd order)
	if (timestepping_order == -2 || timestepping_order == -22)
	{

		////////////////////////////////////////////////////////
		//Calculate order 1 SL-ETDRK (as above) : {U}_1^{n+1} //
		////////////////////////////////////////////////////////
		// Save SL-ETD1RK from above
		SphereData_Spectral U_phi1(sphereDataConfig, 0);
		SphereData_Spectral U_vrt1(sphereDataConfig, 0);
		SphereData_Spectral U_div1(sphereDataConfig, 0);
		U_phi1 = U_phi;
		U_vrt1 = U_vrt;
		U_div1 = U_div;

		/////////////////////////////
		// Apply psi_1 to N(Un) //
		/////////////////////////////
		// already computed in the previous if loop


		//////////////////////////////////////////////////////////////////////////
		// Add half of this to (original) U and interpolate to departure points //
		//////////////////////////////////////////////////////////////////////////
		U_phi = io_U_phi + .5 * i_dt * psi1_FUn_phi;
		U_vrt = io_U_vrt + .5 * i_dt * psi1_FUn_vrt;
		U_div = io_U_div + .5 * i_dt * psi1_FUn_div;

		/////////////////////////////////////////////////////////////
		//Interpolate (.5 * psi1NUn + U) to departure points ( )_* //
		/////////////////////////////////////////////////////////////
		//////SphereData_Spectral psi1_FUn_phi_D(sphereDataConfig, 0);
		//////SphereData_Spectral psi1_FUn_vrt_D(sphereDataConfig, 0);
		//////SphereData_Spectral psi1_FUn_div_D(sphereDataConfig, 0);
		semiLagrangian.apply_sl_timeintegration_uv(
				ops,
				U_phi, U_vrt, U_div,
				pos_lon_d, pos_lat_d,
				psi1_FUn_phi_D, psi1_FUn_vrt_D, psi1_FUn_div_D
			);

		////////////////////////
		// Compute N(U_{n+1}) //
		////////////////////////
		SphereData_Spectral FU1_phi(sphereDataConfig, 0);
		SphereData_Spectral FU1_vrt(sphereDataConfig, 0);
		SphereData_Spectral FU1_div(sphereDataConfig, 0);
		ts_ln_erk_split_uv.euler_timestep_update_lc(
				U_phi1, U_vrt1, U_div1,
				FU1_phi, FU1_vrt, FU1_div,
				i_simulation_timestamp
		);

		if (nonlinear_remainder_treatment == NLRemainderTreatment_enum::NL_REMAINDER_NONLINEAR)
		{
			ts_ln_erk_split_uv.euler_timestep_update_nr(
					U_phi1, U_vrt1, U_div1,
					FU1_phi, FU1_vrt, FU1_div,
					i_simulation_timestamp
			);
		}

		//////////////////////////////
		// Apply psi1 to N(U_{n+1}) //
		//////////////////////////////
		SphereData_Spectral psi1_FU1_phi(sphereDataConfig, 0);
		SphereData_Spectral psi1_FU1_vrt(sphereDataConfig, 0);
		SphereData_Spectral psi1_FU1_div(sphereDataConfig, 0);
		ts_psi1_exp.run_timestep(
				FU1_phi, FU1_vrt, FU1_div,
				psi1_FU1_phi, psi1_FU1_vrt, psi1_FU1_div,
				i_dt,
				i_simulation_timestamp
			);

		/////////////////////////////////////////////////////
		// Add half of this to already interpolated values //
		/////////////////////////////////////////////////////
		U_phi = psi1_FUn_phi_D + .5 * i_dt * psi1_FU1_phi;
		U_vrt = psi1_FUn_vrt_D + .5 * i_dt * psi1_FU1_vrt;
		U_div = psi1_FUn_div_D + .5 * i_dt * psi1_FU1_div;

		////////////////////////////////////////////////////////////////////////
		// Calculate phi_0 of [interpolated (.5 * psi1NUN + U) + .5 * ps1NU1] //
		////////////////////////////////////////////////////////////////////////
		SphereData_Spectral phi0_Un_phi(sphereDataConfig, 0);
		SphereData_Spectral phi0_Un_vrt(sphereDataConfig, 0);
		SphereData_Spectral phi0_Un_div(sphereDataConfig, 0);
		ts_phi0_exp.run_timestep(
				U_phi, U_vrt, U_div,
				phi0_Un_phi, phi0_Un_vrt, phi0_Un_div,
				i_dt,
				i_simulation_timestamp
		);

		U_phi = phi0_Un_phi;
		U_vrt = phi0_Un_vrt;
		U_div = phi0_Un_div;

	}
	if (timestepping_order == -22)
	{

		////////////////////////////////////////////////////////
		//Calculate order 1 SL-ETDRK (as above) : {U}_1^{n+1} //
		////////////////////////////////////////////////////////
		// Save SL-ETD1RK-bis from above
		SphereData_Spectral A_phi1(sphereDataConfig, 0);
		SphereData_Spectral A_vrt1(sphereDataConfig, 0);
		SphereData_Spectral A_div1(sphereDataConfig, 0);
		A_phi1 = U_phi;
		A_vrt1 = U_vrt;
		A_div1 = U_div;

		/////////////////////////////
		// Apply psi_1 to N(Un) //
		/////////////////////////////
		// already computed in the previous if loop

		/////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////// Add half of this to (original) U and interpolate to departure points //
		/////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////U_phi = io_U_phi + .5 * i_dt * psi1_FUn_phi;
		///////////////////////U_vrt = io_U_vrt + .5 * i_dt * psi1_FUn_vrt;
		///////////////////////U_div = io_U_div + .5 * i_dt * psi1_FUn_div;

		////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////Interpolate (.5 * psi1NUn + U) to departure points ( )_* //
		////////////////////////////////////////////////////////////////////////////////////
		///////////////////////SphereData_Spectral psi1_FUn_phi_D(sphereDataConfig, 0);
		///////////////////////SphereData_Spectral psi1_FUn_vrt_D(sphereDataConfig, 0);
		///////////////////////SphereData_Spectral psi1_FUn_div_D(sphereDataConfig, 0);
		///////////////////////semiLagrangian.apply_sl_timeintegration_uv(
		///////////////////////		ops,
		///////////////////////		U_phi, U_vrt, U_div,
		///////////////////////		pos_lon_d, pos_lat_d,
		///////////////////////		psi1_FUn_phi_D, psi1_FUn_vrt_D, psi1_FUn_div_D
		///////////////////////	);

		////////////////////////
		// Compute N(A_{n+1}) //
		////////////////////////
		SphereData_Spectral FA1_phi(sphereDataConfig, 0);
		SphereData_Spectral FA1_vrt(sphereDataConfig, 0);
		SphereData_Spectral FA1_div(sphereDataConfig, 0);
		ts_ln_erk_split_uv.euler_timestep_update_lc(
				A_phi1, A_vrt1, A_div1,
				FA1_phi, FA1_vrt, FA1_div,
				i_simulation_timestamp
		);

		if (nonlinear_remainder_treatment == NLRemainderTreatment_enum::NL_REMAINDER_NONLINEAR)
		{
			ts_ln_erk_split_uv.euler_timestep_update_nr(
					A_phi1, A_vrt1, A_div1,
					FA1_phi, FA1_vrt, FA1_div,
					i_simulation_timestamp
			);
		}

		//////////////////////////////
		// Apply psi1 to N(A_{n+1}) //
		//////////////////////////////
		SphereData_Spectral psi1_FA1_phi(sphereDataConfig, 0);
		SphereData_Spectral psi1_FA1_vrt(sphereDataConfig, 0);
		SphereData_Spectral psi1_FA1_div(sphereDataConfig, 0);
		ts_psi1_exp.run_timestep(
				FA1_phi, FA1_vrt, FA1_div,
				psi1_FA1_phi, psi1_FA1_vrt, psi1_FA1_div,
				i_dt,
				i_simulation_timestamp
			);

		/////////////////////////////////////////////////////
		// Add half of this to already interpolated values //
		/////////////////////////////////////////////////////
		U_phi = psi1_FUn_phi_D + .5 * i_dt * psi1_FA1_phi;
		U_vrt = psi1_FUn_vrt_D + .5 * i_dt * psi1_FA1_vrt;
		U_div = psi1_FUn_div_D + .5 * i_dt * psi1_FA1_div;

		////////////////////////////////////////////////////////////////////////
		// Calculate phi_0 of [interpolated (.5 * psi1NUN + U) + .5 * ps1NU1] //
		////////////////////////////////////////////////////////////////////////
		SphereData_Spectral phi0_Un_phi(sphereDataConfig, 0);
		SphereData_Spectral phi0_Un_vrt(sphereDataConfig, 0);
		SphereData_Spectral phi0_Un_div(sphereDataConfig, 0);
		ts_phi0_exp.run_timestep(
				U_phi, U_vrt, U_div,
				phi0_Un_phi, phi0_Un_vrt, phi0_Un_div,
				i_dt,
				i_simulation_timestamp
		);

		U_phi = phi0_Un_phi;
		U_vrt = phi0_Un_vrt;
		U_div = phi0_Un_div;

	}
	else if (timestepping_order == 2)
	{


		////////////////////////////////////////////////////////
		//Calculate order 1 SL-ETDRK (as above) : {U}_1^{n+1} //
		////////////////////////////////////////////////////////
		// Save SL-ETD1RK from above
		SphereData_Spectral U_phi1(sphereDataConfig, 0);
		SphereData_Spectral U_vrt1(sphereDataConfig, 0);
		SphereData_Spectral U_div1(sphereDataConfig, 0);
		U_phi1 = U_phi;
		U_vrt1 = U_vrt;
		U_div1 = U_div;

		////////////////////////
		// Compute N(U_{n+1}) //
		////////////////////////
		SphereData_Spectral FU1_phi(sphereDataConfig, 0);
		SphereData_Spectral FU1_vrt(sphereDataConfig, 0);
		SphereData_Spectral FU1_div(sphereDataConfig, 0);
		ts_ln_erk_split_uv.euler_timestep_update_lc(
				U_phi1, U_vrt1, U_div1,
				FU1_phi, FU1_vrt, FU1_div,
				i_simulation_timestamp
		);

		if (nonlinear_remainder_treatment == NLRemainderTreatment_enum::NL_REMAINDER_NONLINEAR)
		{
			ts_ln_erk_split_uv.euler_timestep_update_nr(
					U_phi1, U_vrt1, U_div1,
					FU1_phi, FU1_vrt, FU1_div,
					i_simulation_timestamp
			);
		}

		//////////////////////////////
		// Apply psi2 to N(U_{n+1}) //
		//////////////////////////////
		SphereData_Spectral psi2_FU1_phi(sphereDataConfig, 0);
		SphereData_Spectral psi2_FU1_vrt(sphereDataConfig, 0);
		SphereData_Spectral psi2_FU1_div(sphereDataConfig, 0);
		ts_psi2_exp.run_timestep(
				FU1_phi, FU1_vrt, FU1_div,
				psi2_FU1_phi, psi2_FU1_vrt, psi2_FU1_div,
				i_dt,
				i_simulation_timestamp
		);

		//////////////////////////////
		// Apply psi2 to N(U_{n}) //
		//////////////////////////////
		SphereData_Spectral psi2_FUn_phi(sphereDataConfig, 0);
		SphereData_Spectral psi2_FUn_vrt(sphereDataConfig, 0);
		SphereData_Spectral psi2_FUn_div(sphereDataConfig, 0);
		ts_psi2_exp.run_timestep(
				FUn_phi, FUn_vrt, FUn_div,
				psi2_FUn_phi, psi2_FUn_vrt, psi2_FUn_div,
				i_dt,
				i_simulation_timestamp
		);

		////////////////////////////////////////////////////
		// Interpolate psi2NUn to departure points ( )_* //
		////////////////////////////////////////////////////
		SphereData_Spectral psi2_FUn_phi_D(sphereDataConfig, 0);
		SphereData_Spectral psi2_FUn_vrt_D(sphereDataConfig, 0);
		SphereData_Spectral psi2_FUn_div_D(sphereDataConfig, 0);
		semiLagrangian.apply_sl_timeintegration_uv(
				ops,
				psi2_FUn_phi, psi2_FUn_vrt, psi2_FUn_div,
				pos_lon_d, pos_lat_d,
				psi2_FUn_phi_D, psi2_FUn_vrt_D, psi2_FUn_div_D
			);

		////////////////////////////////////////
		// Apply phi0 to psi2NU1 - (psi2NUn)_* //
		////////////////////////////////////////
		SphereData_Spectral phi0_dif2_phi(sphereDataConfig, 0);
		SphereData_Spectral phi0_dif2_vrt(sphereDataConfig, 0);
		SphereData_Spectral phi0_dif2_div(sphereDataConfig, 0);
		ts_phi0_exp.run_timestep(
				psi2_FU1_phi - psi2_FUn_phi_D, psi2_FU1_vrt - psi2_FUn_vrt_D, psi2_FU1_div - psi2_FUn_div_D,
				phi0_dif2_phi, phi0_dif2_vrt, phi0_dif2_div,
				i_dt,
				i_simulation_timestamp
		);

		U_phi = U_phi1 + i_dt*phi0_dif2_phi;
		U_vrt = U_vrt1 + i_dt*phi0_dif2_vrt;
		U_div = U_div1 + i_dt*phi0_dif2_div;


	}
	else if (timestepping_order == 4)
	{

/////////		double dt = i_fixed_dt;
/////////		double dt_half = dt*0.5;
/////////
/////////		/*
/////////		 * Precompute commonly used terms
/////////		 */
/////////		SphereData_Spectral phi0_Un_h(sphereDataConfig);
/////////		SphereData_Spectral phi0_Un_u(sphereDataConfig);
/////////		SphereData_Spectral phi0_Un_v(sphereDataConfig);
/////////
/////////		ts_phi0_rexi.run_timestep(
/////////				io_phi_pert, io_vrt, io_div,
/////////				phi0_Un_h, phi0_Un_u, phi0_Un_v,
/////////				dt_half,
/////////				i_simulation_timestamp
/////////			);
/////////
/////////		SphereData_Spectral FUn_h(sphereDataConfig);
/////////		SphereData_Spectral FUn_u(sphereDataConfig);
/////////		SphereData_Spectral FUn_v(sphereDataConfig);
/////////
/////////		ts_lg_erk_lc_n_erk.euler_timestep_update_lc_n(
/////////				io_phi_pert, io_vrt, io_div,
/////////				FUn_h, FUn_u, FUn_v,
/////////				i_simulation_timestamp
/////////		);
/////////
/////////
/////////
/////////		/*
/////////		 * Some commonly shared buffers
/////////		 */
/////////
/////////		SphereData_Spectral phi1_h(sphereDataConfig);
/////////		SphereData_Spectral phi1_u(sphereDataConfig);
/////////		SphereData_Spectral phi1_v(sphereDataConfig);
/////////
/////////
/////////		/*
/////////		 * A_{n} = \psi_{0}(0.5*\Delta tL)U_{n} + \Delta t\psi_{1}(0.5*\Delta tL) F(U_{n})
/////////		 */
/////////		ts_phi1_rexi.run_timestep(
/////////				FUn_h, FUn_u, FUn_v,
/////////				phi1_h, phi1_u, phi1_v,
/////////				dt_half,
/////////				i_simulation_timestamp
/////////			);
/////////
/////////		SphereData_Spectral A_h = phi0_Un_h + dt_half*phi1_h;
/////////		SphereData_Spectral A_u = phi0_Un_u + dt_half*phi1_u;
/////////		SphereData_Spectral A_v = phi0_Un_v + dt_half*phi1_v;
/////////
/////////
/////////
/////////		/*
/////////		 * B_{n} = \psi_{0}(0.5*\Delta tL)U_{n} + 0.5*\Delta t\psi_{1}(0.5*\Delta tL) F(A_{n}, t_{n} + 0.5*\Delta t)
/////////		 */
/////////
/////////		SphereData_Spectral FAn_h(sphereDataConfig);
/////////		SphereData_Spectral FAn_u(sphereDataConfig);
/////////		SphereData_Spectral FAn_v(sphereDataConfig);
/////////
/////////		ts_lg_erk_lc_n_erk.euler_timestep_update_lc_n(
/////////				A_h, A_u, A_v,
/////////				FAn_h, FAn_u, FAn_v,
/////////				i_simulation_timestamp + dt_half
/////////		);
/////////
/////////		ts_phi1_rexi.run_timestep(
/////////				FAn_h, FAn_u, FAn_v,
/////////				phi1_h, phi1_u, phi1_v,
/////////				dt_half,
/////////				i_simulation_timestamp
/////////			);
/////////
/////////		SphereData_Spectral B_h = phi0_Un_h + dt_half*phi1_h;
/////////		SphereData_Spectral B_u = phi0_Un_u + dt_half*phi1_u;
/////////		SphereData_Spectral B_v = phi0_Un_v + dt_half*phi1_v;
/////////
/////////
/////////
/////////		/*
/////////		 * C_{n} = \psi_{0}(0.5*\Delta tL)U_{n} + 0.5*\Delta t\psi_{1}(0.5* \Delta tL) ( 2 F(B_{n},t_{n} + 0.5*\Delta t)-F(U_{n},t_{n})).
/////////		 */
/////////
/////////		SphereData_Spectral phi0_An_h(sphereDataConfig);
/////////		SphereData_Spectral phi0_An_u(sphereDataConfig);
/////////		SphereData_Spectral phi0_An_v(sphereDataConfig);
/////////
/////////		ts_phi0_rexi.run_timestep(
/////////				A_h, A_u, A_v,
/////////				phi0_An_h, phi0_An_u, phi0_An_v,
/////////				dt_half,
/////////				i_simulation_timestamp
/////////			);
/////////
/////////
/////////		SphereData_Spectral FBn_h(sphereDataConfig);
/////////		SphereData_Spectral FBn_u(sphereDataConfig);
/////////		SphereData_Spectral FBn_v(sphereDataConfig);
/////////
/////////		ts_lg_erk_lc_n_erk.euler_timestep_update_lc_n(
/////////				B_h, B_u, B_v,
/////////				FBn_h, FBn_u, FBn_v,
/////////				i_simulation_timestamp + dt_half
/////////		);
/////////
/////////		ts_phi1_rexi.run_timestep(
/////////				2.0*FBn_h - FUn_h,
/////////				2.0*FBn_u - FUn_u,
/////////				2.0*FBn_v - FUn_v,
/////////				phi1_h,	phi1_u,	phi1_v,
/////////				dt_half,
/////////				i_simulation_timestamp
/////////			);
/////////
/////////		SphereData_Spectral C_h = phi0_An_h + dt_half*phi1_h;
/////////		SphereData_Spectral C_u = phi0_An_u + dt_half*phi1_u;
/////////		SphereData_Spectral C_v = phi0_An_v + dt_half*phi1_v;
/////////
/////////
/////////
/////////		/*
/////////		 * R0 - R3
/////////		 */
/////////		SphereData_Spectral FCn_h(sphereDataConfig);
/////////		SphereData_Spectral FCn_u(sphereDataConfig);
/////////		SphereData_Spectral FCn_v(sphereDataConfig);
/////////
/////////		ts_lg_erk_lc_n_erk.euler_timestep_update_lc_n(
/////////				C_h, C_u, C_v,
/////////				FCn_h, FCn_u, FCn_v,
/////////				i_simulation_timestamp + dt
/////////		);
/////////
/////////		SphereData_Spectral R0_h = io_phi_pert;
/////////		SphereData_Spectral R0_u = io_vrt;
/////////		SphereData_Spectral R0_v = io_div;
/////////
/////////		SphereData_Spectral &R1_h = FUn_h;
/////////		SphereData_Spectral &R1_u = FUn_u;
/////////		SphereData_Spectral &R1_v = FUn_v;
/////////
/////////		SphereData_Spectral R2_h = FAn_h + FBn_h;
/////////		SphereData_Spectral R2_u = FAn_u + FBn_u;
/////////		SphereData_Spectral R2_v = FAn_v + FBn_v;
/////////
/////////		SphereData_Spectral &R3_h = FCn_h;
/////////		SphereData_Spectral &R3_u = FCn_u;
/////////		SphereData_Spectral &R3_v = FCn_v;
/////////
/////////
/////////		/*
/////////		 * U_{n+1} =
/////////		 * 		\psi_{0}(\Delta tL)R_{0}
/////////		 * 			+ \Delta t
/////////		 * 			(
/////////		 * 				  \upsilon_{1}(\Delta tL) R_{1} +
/////////		 * 				2*\upsilon_{2}(\Delta tL) R_{2} +
/////////		 * 				  \upsilon_{3}(\Delta tL) R_{3}
/////////		 * 			)
/////////		 */
/////////		ts_ups0_rexi.run_timestep(
/////////				R0_h, R0_u, R0_v,
/////////				dt,		i_simulation_timestamp
/////////			);
/////////
/////////		ts_ups1_rexi.run_timestep(
/////////				R1_h, R1_u, R1_v,
/////////				dt,		i_simulation_timestamp
/////////			);
/////////
/////////		ts_ups2_rexi.run_timestep(
/////////				R2_h, R2_u, R2_v,
/////////				dt,		i_simulation_timestamp
/////////			);
/////////
/////////		ts_ups3_rexi.run_timestep(
/////////				R3_h, R3_u, R3_v,
/////////				dt,		i_simulation_timestamp
/////////			);
/////////
/////////		io_phi_pert = R0_h + dt*(R1_h + 2.0*R2_h + R3_h);
/////////		io_vrt = R0_u + dt*(R1_u + 2.0*R2_u + R3_u);
/////////		io_div = R0_v + dt*(R1_v + 2.0*R2_v + R3_v);




	}



	// Save current time step for next step
	U_phi_prev = io_U_phi;
	U_vrt_prev = io_U_vrt;
	U_div_prev = io_U_div;

	io_U_phi = U_phi;
	io_U_vrt = U_vrt;
	io_U_div = U_div;

}



/*
 * Setup
 */
void SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etdrk_uv::setup(
		EXP_SimulationVariables &i_rexiSimVars,
		int i_timestepping_order,
		int i_timestepping_order2,
		double i_timestep_size,

		NLRemainderTreatment_enum i_nonlinear_remainder_treatment
)
{
	timestepping_order = i_timestepping_order;
	timestepping_order2 = i_timestepping_order2;

	nonlinear_remainder_treatment = i_nonlinear_remainder_treatment;

	ts_ln_erk_split_uv.setup(i_timestepping_order, true, true, true, true, false);

	// Setup semi-lag
	semiLagrangian.setup(ops.sphereDataConfig);

	if (timestepping_order != timestepping_order2)
		SWEETError("Mismatch of orders, should be equal");

	if (timestepping_order == 0 || timestepping_order == 1 || timestepping_order == -2 || timestepping_order == -22)
	{
		ts_phi0_exp.setup(i_rexiSimVars, "phi0", i_timestep_size, false, true);	/* NO Coriolis */
		ts_psi1_exp.setup(i_rexiSimVars, "psi1", i_timestep_size, false, true);
		ts_chi1_exp.setup(i_rexiSimVars, "chi1", i_timestep_size, false, true, nb_subintegrals);
	}
	else if (timestepping_order == 2)
	{
		ts_phi0_exp.setup(i_rexiSimVars, "phi0", i_timestep_size, false, true);	/* NO Coriolis */
		ts_psi1_exp.setup(i_rexiSimVars, "psi1", i_timestep_size, false, true);
		ts_psi2_exp.setup(i_rexiSimVars, "psi2", i_timestep_size, false, true);
	}
	else if  (timestepping_order == 4)
	{
		SWEETError("4th order method not (yet) supported");

#if 0
		ts_phi0_exp.setup(i_rexiSimVars, "phi0", i_timestep_size*0.5, false, true);	/* NO Coriolis */
		ts_phi1_exp.setup(i_rexiSimVars, "phi1", i_timestep_size*0.5, false, true);
		ts_phi2_exp.setup(i_rexiSimVars, "phi2", i_timestep_size*0.5, false, true);

		// phi0, but with a full time step size
		ts_ups0_exp.setup(i_rexiSimVars, "phi0", i_timestep_size, false, true);
		ts_ups1_exp.setup(i_rexiSimVars, "ups1", i_timestep_size, false, true);
		ts_ups2_exp.setup(i_rexiSimVars, "ups2", i_timestep_size, false, true);
		ts_ups3_exp.setup(i_rexiSimVars, "ups3", i_timestep_size, false, true);
#endif

	}
	else
	{
		SWEETError("TODO: This order is not implemented, yet!");
	}
}


SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etdrk_uv::SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etdrk_uv(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		ops(i_op),

		ts_ln_erk_split_uv(simVars, ops),

		ts_phi0_exp(simVars, ops),
		ts_phi2_exp(simVars, ops),

		ts_psi1_exp(simVars, ops),
		ts_psi2_exp(simVars, ops),

		ts_chi1_exp(simVars, ops),

		semiLagrangian(simVars),
		sphereSampler(semiLagrangian.sphereSampler)

#if 0
		,
		ts_ups0_exp(simVars, ops),
		ts_ups1_exp(simVars, ops),
		ts_ups2_exp(simVars, ops),
		ts_ups3_exp(simVars, ops)
#endif
{
}



SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etdrk_uv::~SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etdrk_uv()
{
}

