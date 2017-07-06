/*
 * PlaneOperatorsComplex.hpp
 *
 *  Created on: 30 Jun 2015
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */
#ifndef SRC_INCLUDE_SWEET_PLANE_OPERATORS_COMPLEX_HPP_
#define SRC_INCLUDE_SWEET_PLANE_OPERATORS_COMPLEX_HPP_


#include <sweet/plane/PlaneDataComplex.hpp>
#include <sweet/plane/PlaneDataConfig.hpp>


class PlaneOperatorsComplex
{
	PlaneDataConfig *planeDataConfig;

public:
	// differential operators
	PlaneDataComplex diff_c_x, diff_c_y;
	PlaneDataComplex diff2_c_x, diff2_c_y;

	/**
	 * D2, e.g. for viscosity
	 */
	PlaneDataComplex diff2(
			const PlaneDataComplex &i_dataArray
	)
	{
		return diff2_c_x(i_dataArray) + diff2_c_y(i_dataArray);
	}



	/**
	 *        __2
	 * apply  \/  operator (aka Laplace)
	 */
	inline PlaneDataComplex laplace(
			const PlaneDataComplex &i_a
	)
	{
		return diff2_c_x(i_a)+diff2_c_y(i_a);
	}


	/**
	 *        __
	 * apply  \/ .  operator
	 */
	inline PlaneDataComplex diff_dot(
			const PlaneDataComplex &i_a
	)
	{
		return diff_c_x(i_a)+diff_c_y(i_a);
	}



	/**
	 * Diff N operator for hyperviscosity, see
	 * "Numerical Techniques for Global Atmospheric Models", page 500
	 */
	inline PlaneDataComplex diffN_x(
			const PlaneDataComplex &io_u,
			int i_order
	)
	{
		if (i_order == 0)
			return io_u;

		PlaneDataComplex tu = io_u;

		for (int i = 0; i < i_order/2; i++)
			tu = diff2_c_x(tu);

		if (i_order & 1)
			tu = diff_c_x(tu);

		return tu;
	}


	/**
	 * Diff N operator for hyperviscosity, see
	 * "Numerical Techniques for Global Atmospheric Models", page 500
	 */
	inline PlaneDataComplex diffN_y(
			const PlaneDataComplex &io_v,
			int i_order
	)
	{
		if (i_order == 0)
			return io_v;

		PlaneDataComplex tv = io_v;

		for (int i = 0; i < i_order/2; i++)
			tv = diff2_c_y(tv);

		if (i_order & 1)
			tv = diff_c_y(tv);

		return tv;
	}


#if SWEET_USE_PLANE_SPECTRAL_SPACE
	/**
	 * Diffusion or hyperviscosity coefficients
	 * Simply calculates the spectral coefficients
	 *  see "Numerical Techniques for Global Atmospheric Models", page 500
	 *
	 * i_order (q) needs to be even!!! (second or forth order usually)
	 *
	 * Returns operator D^q
	 *
	 */
	inline PlaneDataComplex diffusion_coefficient(
			int i_order
	)
	{
		//Check if even
		assert( i_order % 2 == 0);
		assert( i_order > 0);
		PlaneDataComplex out = diff2_c_x+diff2_c_y;

		for (int i = 1; i < i_order/2; i++)
			out = pow(-1, i)*(diff2_c_x(out)+diff2_c_y(out));

		return out;
	}

	/**
	 * Calculates implicit diffusion (applies 1/(1-mu*dt*D^q) to spectrum)
	 *  see "Numerical Techniques for Global Atmospheric Models", page 500
	 *
	 * i_order (q) needs to be even!!! (second or forth order usually)
	 * i_coef is mu*dt
	 *
	 * Only works in spectral space
	 *
	 */
	inline PlaneDataComplex implicit_diffusion(
			const PlaneDataComplex &i_data,
			double i_coef,
			int i_order
	)
	{
		PlaneDataComplex out=i_data;

		// Get diffusion coefficients (these are the -mu*dt*D^q, where q is the order
		PlaneDataComplex diff = -i_coef*diffusion_coefficient(i_order);

		// Add 1 to get denominator
		diff = diff.spectral_addScalarAll(1.0);

		// Invert
		diff = diff.spectral_invert();
		// apply to data
		out=diff(out);
		return out;
	}
#endif


	void setup(
			PlaneDataConfig *i_planeDataConfig,
			const double i_domain_size[2]
	)
	{
		planeDataConfig = i_planeDataConfig;
		
		diff_c_x.setup(i_planeDataConfig);
		diff_c_y.setup(i_planeDataConfig);
		diff2_c_x.setup(i_planeDataConfig);
		diff2_c_y.setup(i_planeDataConfig);
		

		/*
		 * setup spectral differential operators
		 * 		diff(e(ix), x)
		 */
		// Assume, that errors are linearly depending on the resolution
		// see test_spectral_ops.cpp

#if !SWEET_USE_PLANE_SPECTRAL_SPACE
		std::cerr << "Activate spectral space during compile time to use spectral diffs. Otherwise, the convolution would be freakingly expensive" << std::endl;
		assert(false);
#endif

		{
		        diff_c_x.spectral_set_all(0, 0);
			double scale = 2.0*M_PIl/i_domain_size[0];

//#if !SWEET_REXI_THREAD_PARALLEL_SUM
//#		pragma omp parallel for OPENMP_PAR_SIMD
//#endif
			for (std::size_t j = 0; j < planeDataConfig->spectral_complex_data_size[1]/2; j++)
			{
				for (std::size_t i = 1; i < planeDataConfig->spectral_complex_data_size[0]/2; i++)
				{
					diff_c_x.p_spectral_set(
							j,
							i,
							0,
							(double)i*scale
						);
					diff_c_x.p_spectral_set(
							planeDataConfig->spectral_complex_data_size[1]-1-j,
							i,
							0,
							(double)i*scale
						);

					diff_c_x.p_spectral_set(
							j,
							planeDataConfig->spectral_complex_data_size[0]-i,
							0,
							-(double)i*scale
						);
					diff_c_x.p_spectral_set(
							planeDataConfig->spectral_complex_data_size[1]-1-j,
							planeDataConfig->spectral_complex_data_size[0]-i,
							0,
							-(double)i*scale
						);
				}
			}

			diff_c_x.physical_space_data_valid = false;
			diff_c_x.spectral_space_data_valid = true;
		}


		/*
		 * DIFF operator in y axis
		 */
		{
			diff_c_y.spectral_set_all(0, 0);

			double scale = 2.0*M_PIl/i_domain_size[1];

//#if !SWEET_REXI_THREAD_PARALLEL_SUM
//#		pragma omp parallel for OPENMP_PAR_SIMD
//#endif
			for (int j = 1; j < (int)planeDataConfig->spectral_complex_data_size[1]/2; j++)
			{
				for (int i = 0; i < (int)planeDataConfig->spectral_complex_data_size[0]/2; i++)
				{
					diff_c_y.p_spectral_set(
							j,
							i,
							0,
							(double)j*scale
						);
					diff_c_y.p_spectral_set(
							planeDataConfig->spectral_complex_data_size[1]-j,
							i,
							0,
							-(double)j*scale
						);

					diff_c_y.p_spectral_set(
							j,
							planeDataConfig->spectral_complex_data_size[0]-i-1,
							0,
							(double)j*scale
						);
					diff_c_y.p_spectral_set(
							planeDataConfig->spectral_complex_data_size[1]-j,
							planeDataConfig->spectral_complex_data_size[0]-i-1,
							0,
							-(double)j*scale
						);
				}
			}

			diff_c_y.physical_space_data_valid = false;
			diff_c_y.spectral_space_data_valid = true;
		}

		/*
		 * 2nd order differential operators
		 */
		{
			diff2_c_x.spectral_set_all(0, 0);

			double scale = 2.0*M_PIl/i_domain_size[0];

//#if !SWEET_REXI_THREAD_PARALLEL_SUM
//#		pragma omp parallel for OPENMP_PAR_SIMD
//#endif
			for (std::size_t j = 0; j < planeDataConfig->spectral_complex_data_size[1]/2; j++)
			{
				for (std::size_t i = 1; i < planeDataConfig->spectral_complex_data_size[0]/2; i++)
				{
					diff2_c_x.p_spectral_set(
							j,
							i,
							-(double)i*scale*(double)i*scale,
							0
						);
					diff2_c_x.p_spectral_set(
							planeDataConfig->spectral_complex_data_size[1]-1-j,
							i,
							-(double)i*scale*(double)i*scale,
							0
						);

					diff2_c_x.p_spectral_set(	j,
							planeDataConfig->spectral_complex_data_size[0]-i,
							-(double)i*scale*(double)i*scale,
							0
						);
					diff2_c_x.p_spectral_set(
							planeDataConfig->spectral_complex_data_size[1]-1-j,
							planeDataConfig->spectral_complex_data_size[0]-i,
							-(double)i*scale*(double)i*scale,
							0
						);
				}
			}

			diff2_c_x.physical_space_data_valid = false;
			diff2_c_x.spectral_space_data_valid = true;
		}

		{
			diff2_c_y.spectral_set_all(0, 0);
			double scale = 2.0*M_PIl/i_domain_size[1];

			//#if !SWEET_REXI_THREAD_PARALLEL_SUM
			//#		pragma omp parallel for OPENMP_PAR_SIMD
			//#endif
			for (int j = 1; j < (int)planeDataConfig->spectral_complex_data_size[1]/2; j++)
			{
				for (int i = 0; i < (int)planeDataConfig->spectral_complex_data_size[0]/2; i++)
				{
					diff2_c_y.p_spectral_set(
							j,
							i,
							-(double)j*scale*(double)j*scale,
							0
						);
					diff2_c_y.p_spectral_set(
							planeDataConfig->spectral_complex_data_size[1]-j,
							i,
							-(double)j*scale*(double)j*scale,
							0
						);

					diff2_c_y.p_spectral_set(
							j,
							planeDataConfig->spectral_complex_data_size[0]-i-1,
							-(double)j*scale*(double)j*scale,
							0
						);
					diff2_c_y.p_spectral_set(
							planeDataConfig->spectral_complex_data_size[1]-j,
							planeDataConfig->spectral_complex_data_size[0]-i-1,
							-(double)j*scale*(double)j*scale,
							0
						);
				}
			}

			diff2_c_y.physical_space_data_valid = false;
			diff2_c_y.spectral_space_data_valid = true;
		}

	}



public:
	PlaneOperatorsComplex()
	{
	}



public:
	PlaneOperatorsComplex(
		PlaneDataConfig *i_planeDataConfig,
		const double i_domain_size[2]	///< domain size
	)	:
		planeDataConfig(i_planeDataConfig),

		diff_c_x(i_planeDataConfig),
		diff_c_y(i_planeDataConfig),
		diff2_c_x(i_planeDataConfig),
		diff2_c_y(i_planeDataConfig)
	{
		setup(i_planeDataConfig, i_domain_size);
	}
};



#endif
