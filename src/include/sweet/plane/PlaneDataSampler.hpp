/*
 * Sampler2D.hpp
 *
 *  Created on: 4 Dec 2015
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */
#ifndef SRC_INCLUDE_SWEET_PLANEDATASAMPLER_HPP_
#define SRC_INCLUDE_SWEET_PLANEDATASAMPLER_HPP_

#include <sweet/ScalarDataArray.hpp>
//#include "PlaneDataComplex.hpp"


/**
 * this is a sampler class which provides method to provide
 * interpolated sampled values on 2D physical data which
 * is provided by PlaneData
 */
class PlaneDataSampler
{
public:
	double domain_size[2];			/// real physical size of the domain
	int res[2];						/// resolution of domain
	const PlaneDataConfig *planeDataConfig;

private:
	double cached_scale_factor[2];			/// cached parameters for sampling


public:
	PlaneDataSampler(
		double i_domain_size[2],	/// real physical size of the domain
		const PlaneDataConfig *i_planeDataConfig
	)
	{
		assert(i_planeDataConfig != nullptr);

		planeDataConfig = i_planeDataConfig;
		setup(i_domain_size, planeDataConfig);
	}


	PlaneDataSampler()
	{
		planeDataConfig = nullptr;

		res[0] = -1;
		res[1] = -1;

		cached_scale_factor[0] = -1;
		cached_scale_factor[1] = -1;

		domain_size[0] = -1;
		domain_size[1] = -1;
	}



public:
	void setup(
		double i_domain_size[2],	/// real physical size of the domain
		const PlaneDataConfig *i_planeDataConfig
	)
	{
		assert(i_planeDataConfig != nullptr);
		planeDataConfig = i_planeDataConfig;

		domain_size[0] = i_domain_size[0];
		domain_size[1] = i_domain_size[1];

		res[0] = i_planeDataConfig->physical_res[0];
		res[1] = i_planeDataConfig->physical_res[1];

		cached_scale_factor[0] = (double)i_planeDataConfig->physical_res[0] / i_domain_size[0];
		cached_scale_factor[1] = (double)i_planeDataConfig->physical_res[1] / i_domain_size[1];
	}

public:
	/**
	 * wrap the position i in a periodic domain of size i_res
	 */
#if 0

#error "This modulo operation doesn't work!"
	inline
	int wrapPeriodic(int i, int i_res)
	{
		return (i + i_res*10) % i_res;
	}

	inline
	double wrapPeriodic(double i, double i_res)
	{
		return fmodf(i + i_res*10.0f, i_res);
	}

#else

	template <typename T>
	inline
	double wrapPeriodic(T i, T i_res)
	{
#if 1
		int c = 10;
		while (i < 0 && c-- > 0)
			i += i_res;

		int d = 10;
		while (i >= i_res && d-- > 0)
			i -= i_res;
#elif 1

			i = (i + i_res*10) % i_res;
		else if (typeid(T) == typeid(double))
			i = fmod(i + i_res*10.0, i_res);
		else
			i = fmodf(i + i_res*10.0f, i_res);

#else
		if (i < 0)
			i += i_res;

		if (i >= i_res)
			i -= i_res;
#endif

		if (i < 0 || i >= i_res)
			SWEETError("Stopping here: Probably an unstable velocity field since more than one periodic movement exists.");

		assert(i >= 0 && i < i_res);

		return i;
	}
#endif


public:
	void bicubic_scalar(
			const PlaneData_Physical  &i_data,			///< sampling data

			const ScalarDataArray &i_pos_x,		///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,		///< y positions of interpolation points

			double *o_data,						///< output values

			double i_shift_x = 0.0,				///< shift in x for staggered grids
			double i_shift_y = 0.0				///< shift in y for staggered grids
	)
	{
		assert(res[0] > 0);
		assert(cached_scale_factor[0] > 0);
		assert(i_pos_x.number_of_elements == i_pos_y.number_of_elements);

		std::size_t max_pos_idx = i_pos_x.number_of_elements;

#if SWEET_DEBUG
#if SWEET_THREADING_SPACE || SWEET_THREADING_TIME_REXI

		if (omp_get_num_threads() > 1)
			SWEETError("Are we already in parallel region? Threading race conditions likely!");
#endif
#endif

		// iterate over all positions in parallel
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t pos_idx = 0; pos_idx < max_pos_idx; pos_idx++)
		{
			/*
			 * load position to interpolate
			 * posx stores all x-coordinates of the arrival points
			 * posy stores all y-coordinates of the arrival points
			 *
			 * Both are arrays and matching array indices (pos_idx) below index the coordinates for the same point.
			 *
			 * Scale factor (Nx/dx, Ny/dy) maps from the physical space to the array space.
			 * The array space is from [0; N[
			 *
			 * shift_x/y is operating in array space. Hence, staggered grid can be
			 * implemented by setting this to 0.5 or -0.5
			 * for C grid, to interpolate given u data, use i_shift_x = 0.0,  i_shift_y = -0.5
			 *             to interpolate given v data, use i_shift_y = -0.5, i_shift_y = 0.0
			 *  pay attention to the negative shift, which is necessary because the staggered grids are positively shifted
			 *  and this shift has to be removed for the interpolation
			 */
			double pos_x = wrapPeriodic(i_pos_x.scalar_data[pos_idx]*cached_scale_factor[0] + i_shift_x, (double)res[0]);
			double pos_y = wrapPeriodic(i_pos_y.scalar_data[pos_idx]*cached_scale_factor[1] + i_shift_y, (double)res[1]);

			/**
			 * For the interpolation, we assume node-aligned values
			 *
			 * x0  x1  x2  x3
			 * |---|---|---|---
			 * 0   2   4   6    <- positions and associated values e.g. for domain size 8
			 */

			/**
			 * See http://www.paulinternet.nl/?page=bicubic
			 */
			// compute x/y position
			double x = pos_x - floor(pos_x);
			double y = pos_y - floor(pos_y);

			// precompute x-position indices since they are reused 4 times
			int idx_i[4];
			{

// TODO: Each 2nd wrapPeriodic is obsolete!!!
				int i = (int)pos_x-1;

				i = wrapPeriodic(i, res[0]);
				idx_i[0] = wrapPeriodic(i, res[0]);

				i = wrapPeriodic(i+1, res[0]);
				idx_i[1] = wrapPeriodic(i, res[0]);

				i = wrapPeriodic(i+1, res[0]);
				idx_i[2] = wrapPeriodic(i, res[0]);

				i = wrapPeriodic(i+1, res[0]);
				idx_i[3] = wrapPeriodic(i, res[0]);
			}

			/**
			 * iterate over rows and interpolate over the columns in the x direction
			 */
			// start at this row
			int idx_j = wrapPeriodic((int)pos_y-1, res[1]);

			double q[4];
			for (int kj = 0; kj < 4; kj++)
			{
				double p[4];

				p[0] = i_data.physical_get(idx_j, idx_i[0]);
				p[1] = i_data.physical_get(idx_j, idx_i[1]);
				p[2] = i_data.physical_get(idx_j, idx_i[2]);
				p[3] = i_data.physical_get(idx_j, idx_i[3]);

				q[kj] = p[1] + 0.5 * x*(p[2] - p[0] + x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0])));

				idx_j = wrapPeriodic(idx_j+1, res[1]);
			}
			double value = q[1] + 0.5 * y*(q[2] - q[0] + y*(2.0*q[0] - 5.0*q[1] + 4.0*q[2] - q[3] + y*(3.0*(q[1] - q[2]) + q[3] - q[0])));

			o_data[pos_idx] = value;
		}
	}



public:
	void bicubic_scalar(
			const PlaneData_Physical  &i_data,			///< sampling data

			const ScalarDataArray &i_pos_x,		///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,		///< y positions of interpolation points

			PlaneData_Physical  &o_data,					///< output values

			double i_shift_x = 0.0,				///< shift in x for staggered grids
			double i_shift_y = 0.0				///< shift in y for staggered grids
	)
	{
		assert(res[0] > 0);
		assert(cached_scale_factor[0] > 0);
		assert(i_pos_x.number_of_elements == i_pos_y.number_of_elements);
		assert(i_pos_x.number_of_elements == o_data.planeDataConfig->physical_array_data_number_of_elements);

		bicubic_scalar(i_data, i_pos_x, i_pos_y, o_data.physical_space_data, i_shift_x, i_shift_y);


	}


public:
	void bicubic_scalar(
			const PlaneData_Physical  &i_data,				///< sampling data

			const ScalarDataArray &i_pos_x,			///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,			///< y positions of interpolation points

			ScalarDataArray &o_data,				///< output values

			double i_shift_x = 0.0,            ///< shift in x for staggered grids
			double i_shift_y = 0.0				///< shift in y for staggered grids
	)
	{
		assert(res[0] > 0);
		assert(cached_scale_factor[0] > 0);
		assert(i_pos_x.number_of_elements == i_pos_y.number_of_elements);
		assert(i_pos_x.number_of_elements == o_data.number_of_elements);

		bicubic_scalar(i_data, i_pos_x, i_pos_y, o_data.scalar_data, i_shift_x, i_shift_y);
	}


public:
	void bilinear_scalar(
			const PlaneData_Physical  &i_data,				///< sampling data

			const ScalarDataArray &i_pos_x,				///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,				///< y positions of interpolation points

			double *o_data,							///< output values
			double i_shift_x = 0.0,
			double i_shift_y = 0.0
	)
	{
		/*
		 * SHIFT - important
		 * for C grid, to interpolate given u data, use i_shift_x = 0.0,  i_shift_y = -0.5
		 *             to interpolate given v data, use i_shift_y = -0.5, i_shift_y = 0.0
		 *  pay attention to the negative shift, which is necessary because the staggered grids are positively shifted
		 *  and this shift has to be removed for the interpolation
		 */


		std::size_t size = i_pos_x.number_of_elements;

		assert(size != 0);

		// iterate over all positions
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t pos_idx = 0; pos_idx < size; pos_idx++)
		{
			// load position to interpolate
			double pos_x = wrapPeriodic(i_pos_x.scalar_data[pos_idx]*cached_scale_factor[0] + i_shift_x, (double)res[0]);
			double pos_y = wrapPeriodic(i_pos_y.scalar_data[pos_idx]*cached_scale_factor[1] + i_shift_y, (double)res[1]);

			/**
			 * See http://www.paulinternet.nl/?page=bicubic
			 */
			// compute x/y position
			double x = pos_x - floor(pos_x);
			double y = pos_y - floor(pos_y);

			// precompute x-position indices since they are reused 2 times
			int idx_i[2];
			{
				int i = (int)pos_x;
				idx_i[0] = i;

// TODO: Each 2nd wrapPeriodic is obsolete!!!
				i = wrapPeriodic(i+1, res[0]);
				idx_i[1] = wrapPeriodic(i, res[0]);
			}

			/**
			 * iterate over rows and interpolate over the columns in the x direction
			 */
			// start at this row
			int idx_j = pos_y;//wrapPeriodic((int)pos_y, res[1]);

			double q[2];
			for (int kj = 0; kj < 2; kj++)
			{
				double p[2];
				p[0] = i_data.physical_get(idx_j, idx_i[0]);
				p[1] = i_data.physical_get(idx_j, idx_i[1]);

				q[kj] = p[0] + x*(p[1]-p[0]);

				idx_j = wrapPeriodic(idx_j+1, res[1]);
				//std::cout<< p[0] << p[1] << x << std::endl;
			}

			double value = q[0] + y*(q[1]-q[0]);

			o_data[pos_idx] = value;
		}
	}


public:
	void bilinear_scalar(
			const PlaneData_Physical  &i_data,				///< sampling data

			const ScalarDataArray &i_pos_x,				///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,				///< y positions of interpolation points

			ScalarDataArray &o_data,				///< output values
			double i_shift_x = 0.0,
			double i_shift_y = 0.0
	)
	{
		bilinear_scalar(i_data, i_pos_x, i_pos_y, o_data.scalar_data, i_shift_x, i_shift_y);
	}


public:
	void bilinear_scalar(
			const PlaneData_Physical  &i_data,				///< sampling data

			const ScalarDataArray &i_pos_x,				///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,				///< y positions of interpolation points

			PlaneData_Physical  &o_data,				///< output values

			double i_shift_x = 0.0,
			double i_shift_y = 0.0
	)
	{
		/*
		 * SHIFT - important
		 * for C grid, to interpolate given u data, use i_shift_x = 0.0,  i_shift_y = -0.5
		 *             to interpolate given v data, use i_shift_y = -0.5, i_shift_y = 0.0
		 *  pay attention to the negative shift, which is necessary because the staggered grids are positively shifted
		 *  and this shift has to be removed for the interpolation
		 */
		assert(i_pos_x.number_of_elements == i_pos_y.number_of_elements);
		assert(i_pos_x.number_of_elements == o_data.planeDataConfig->physical_array_data_number_of_elements);

		bilinear_scalar(i_data, i_pos_x, i_pos_y, o_data.physical_space_data, i_shift_x, i_shift_y);

	}


public:
	const ScalarDataArray bilinear_scalar(
			const PlaneData_Physical  &i_data,				///< sampling data

			const ScalarDataArray &i_pos_x,				///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,				///< y positions of interpolation points

			double i_shift_x = 0.0,
			double i_shift_y = 0.0
	)
	{
		ScalarDataArray out(i_data.planeDataConfig->physical_array_data_number_of_elements);
		bilinear_scalar(i_data, i_pos_x, i_pos_y, out, i_shift_x, i_shift_y);
		return out;
	}

public:
	const PlaneData_Physical  bicubic_scalar(
			PlaneData_Physical  &i_data,				///< sampling data

			const ScalarDataArray &i_pos_x,				///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,				///< y positions of interpolation points
			//PlaneData* i_pos[2],	///< sampling position
			double i_shift_x = 0.0,
			double i_shift_y = 0.0
	)
	{
		PlaneData_Physical  out(planeDataConfig);
		bicubic_scalar(i_data, i_pos_x, i_pos_y, out, i_shift_x, i_shift_y);
		return out;
	}


// Same interface functions but with PlaneData_Spectral as argument

public:
	void bilinear_scalar(
			const PlaneData_Spectral &i_data,				///< sampling data

			const ScalarDataArray &i_pos_x,				///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,				///< y positions of interpolation points

			ScalarDataArray &o_data,				///< output values
			double i_shift_x = 0.0,
			double i_shift_y = 0.0
	)
	{
		PlaneData_Physical  tmp = i_data.toPhys();
		bilinear_scalar(tmp, i_pos_x, i_pos_y, o_data.scalar_data, i_shift_x, i_shift_y);
	}


public:
	void bilinear_scalar(
			const PlaneData_Spectral &i_data,				///< sampling data

			const ScalarDataArray &i_pos_x,				///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,				///< y positions of interpolation points

			PlaneData_Spectral &o_data,				///< output values

			double i_shift_x = 0.0,
			double i_shift_y = 0.0
	)
	{

		PlaneData_Physical  tmp = i_data.toPhys();
		PlaneData_Physical  o_data_phys(i_data.planeDataConfig);

		/*
		 * SHIFT - important
		 * for C grid, to interpolate given u data, use i_shift_x = 0.0,  i_shift_y = -0.5
		 *             to interpolate given v data, use i_shift_y = -0.5, i_shift_y = 0.0
		 *  pay attention to the negative shift, which is necessary because the staggered grids are positively shifted
		 *  and this shift has to be removed for the interpolation
		 */
		assert(i_pos_x.number_of_elements == i_pos_y.number_of_elements);
		assert(i_pos_x.number_of_elements == o_data.planeDataConfig->physical_array_data_number_of_elements);

		bilinear_scalar(tmp, i_pos_x, i_pos_y, o_data_phys.physical_space_data, i_shift_x, i_shift_y);

		o_data.loadPlaneDataPhysical(o_data_phys);

	}




public:
	const PlaneData_Physical  bicubic_scalar(
			PlaneData_Spectral &i_data,				///< sampling data

			const ScalarDataArray &i_pos_x,				///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,				///< y positions of interpolation points
			//PlaneData* i_pos[2],	///< sampling position
			double i_shift_x = 0.0,
			double i_shift_y = 0.0
	)
	{
		PlaneData_Physical  tmp = i_data.toPhys();
		PlaneData_Physical  out(i_data.planeDataConfig);
		bicubic_scalar(tmp, i_pos_x, i_pos_y, out, i_shift_x, i_shift_y);
		return out;
	}

public:
	void bicubic_scalar(
			const PlaneData_Spectral &i_data,			///< sampling data

			const ScalarDataArray &i_pos_x,		///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,		///< y positions of interpolation points

			PlaneData_Spectral &o_data,					///< output values

			double i_shift_x = 0.0,				///< shift in x for staggered grids
			double i_shift_y = 0.0				///< shift in y for staggered grids
	)
	{
		PlaneData_Physical  tmp = i_data.toPhys();
		PlaneData_Physical  o_data_phys(i_data.planeDataConfig);

		assert(res[0] > 0);
		assert(cached_scale_factor[0] > 0);
		assert(i_pos_x.number_of_elements == i_pos_y.number_of_elements);
		assert(i_pos_x.number_of_elements == o_data.planeDataConfig->physical_array_data_number_of_elements);

		bicubic_scalar(tmp, i_pos_x, i_pos_y, o_data_phys.physical_space_data, i_shift_x, i_shift_y);

		o_data.loadPlaneDataPhysical(o_data_phys);

	}



public:
	void bicubic_scalar_NEW2(
			const PlaneData_Physical &i_data,			///< sampling data

			const ScalarDataArray &i_pos_x,		///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,		///< y positions of interpolation points

			double *o_data,						///< output values

			double i_shift_x = 0.0,				///< shift in x for staggered grids
			double i_shift_y = 0.0				///< shift in y for staggered grids
	)
	{
		assert(res[0] > 0);
		assert(cached_scale_factor[0] > 0);
		assert(i_pos_x.number_of_elements == i_pos_y.number_of_elements);

		std::size_t max_pos_idx = i_pos_x.number_of_elements;

#if SWEET_DEBUG
#if SWEET_THREADING_SPACE || SWEET_THREADING_TIME_REXI

		if (omp_get_num_threads() > 1)
			SWEETError("Are we already in parallel region? Threading race conditions likely!");
#endif
#endif

		// iterate over all positions in parallel
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t pos_idx = 0; pos_idx < max_pos_idx; pos_idx++)
		{
			/*
			 * load position to interpolate
			 * posx stores all x-coordinates of the arrival points
			 * posy stores all y-coordinates of the arrival points
			 *
			 * Both are arrays and matching array indices (pos_idx) below index the coordinates for the same point.
			 
			 * Scale factor (Nx/dx, Ny/dy) maps from the physical space to the array space.
			 * The array space is from [0; N[
			 *
			 * shift_x/y is operating in array space. Hence, staggered grid can be
			 * implemented by setting this to 0.5 or -0.5
			 * for C grid, to interpolate given u data, use i_shift_x = 0.0,  i_shift_y = -0.5
			 *             to interpolate given v data, use i_shift_y = -0.5, i_shift_y = 0.0
			 *  pay attention to the negative shift, which is necessary because the staggered grids are positively shifted
			 *  and this shift has to be removed for the interpolation
			 */
			double pos_x = wrapPeriodic(i_pos_x.scalar_data[pos_idx]*cached_scale_factor[0] + i_shift_x, (double)res[0]);
			double pos_y = wrapPeriodic(i_pos_y.scalar_data[pos_idx]*cached_scale_factor[1] + i_shift_y, (double)res[1]);

			/**
			 * For the interpolation, we assume node-aligned values
			 *
			 * x0  x1  x2  x3
			 * |---|---|---|---
			 * 0   2   4   6    <- positions and associated values e.g. for domain size 8
			 */

			/**
			 * See http://www.paulinternet.nl/?page=bicubic
			 */
			// compute x/y position
			double x = pos_x - floor(pos_x);
			double y = pos_y - floor(pos_y);

			// precompute x-position indices since they are reused 4 times
			int idx_i[10];  // up to 8th order center approximation to first derivative
			{

// TODO: Each 2nd wrapPeriodic is obsolete!!!
				int i = (int)pos_x-4; // C2
				/////////////int i = (int)pos_x-2;

				i = wrapPeriodic(i, res[0]);
				idx_i[0] = wrapPeriodic(i, res[0]);

                                for (int j = 1; j < 10; j++)
				{
					i = wrapPeriodic(i + 1, res[0]);
					idx_i[j] = wrapPeriodic(i, res[0]);
				}

			}

			/**
			 * iterate over rows and interpolate over the columns in the x direction
			 */
			// start at this row
			/////////////////int idx_j = wrapPeriodic((int)pos_y-2, res[1]);
			int idx_j = wrapPeriodic((int)pos_y-4, res[1]);

			double q[10];
			for (int kj = 0; kj < 10; kj++)
			{
				double p[10];

                                for (int j = 0; j < 10; j++)
					p[j] = i_data.physical_get(idx_j, idx_i[j]);

				double x0 = 0.;
				double x1 = 1.;
				double p0 = p[4];
				double p1 = p[5];

				double derx0  = (p[5] - p[3]) / 2.;
				double derx1  = (p[6] - p[4]) / 2.;

				double xx0 = x - x0;
				double xx0_2 = xx0 * xx0;
				double xx1 = x - x1;
				double x1x0 = x1 - x0;
	
				q[kj] = p0 + derx0 * xx0
					   + (p1 - p0 - derx0 * x1x0) * xx0_2 / (x1x0 * x1x0)
					   + ( (derx0 + derx1) * x1x0 - 2. * (p1 - p0) ) * xx0_2 * xx1 / (x1x0 * x1x0 * x1x0);

				idx_j = wrapPeriodic(idx_j+1, res[1]);
			}

			double y0 = 0.;
			double y1 = 1.;
			double q0 = q[4];
			double q1 = q[5];

			double dery0  = (q[5] - q[3]) / 2.;
			double dery1  = (q[6] - q[4]) / 2.;

			double yy0 = y - y0;
			double yy0_2 = yy0 * yy0;
			double yy1 = y - y1;
			double y1y0 = y1 - y0;

			double value = q0 + dery0 * yy0
					  + (q1 - q0 - dery0 * y1y0) * yy0_2 / (y1y0 * y1y0)
					  + ( (dery0 + dery1) * y1y0 - 2. * (q1 - q0) ) * yy0_2 * yy1 / (y1y0 * y1y0 * y1y0);

			o_data[pos_idx] = value;
		}
	}


public:
	void bicubic_scalar_NEW4(
			const PlaneData_Physical &i_data,			///< sampling data

			const ScalarDataArray &i_pos_x,		///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,		///< y positions of interpolation points

			double *o_data,						///< output values

			double i_shift_x = 0.0,				///< shift in x for staggered grids
			double i_shift_y = 0.0				///< shift in y for staggered grids
	)
	{
		assert(res[0] > 0);
		assert(cached_scale_factor[0] > 0);
		assert(i_pos_x.number_of_elements == i_pos_y.number_of_elements);

		std::size_t max_pos_idx = i_pos_x.number_of_elements;

#if SWEET_DEBUG
#if SWEET_THREADING_SPACE || SWEET_THREADING_TIME_REXI

		if (omp_get_num_threads() > 1)
			SWEETError("Are we already in parallel region? Threading race conditions likely!");
#endif
#endif

		// iterate over all positions in parallel
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t pos_idx = 0; pos_idx < max_pos_idx; pos_idx++)
		{
			/*
			 * load position to interpolate
			 * posx stores all x-coordinates of the arrival points
			 * posy stores all y-coordinates of the arrival points
			 *
			 * Both are arrays and matching array indices (pos_idx) below index the coordinates for the same point.
			 
			 * Scale factor (Nx/dx, Ny/dy) maps from the physical space to the array space.
			 * The array space is from [0; N[
			 *
			 * shift_x/y is operating in array space. Hence, staggered grid can be
			 * implemented by setting this to 0.5 or -0.5
			 * for C grid, to interpolate given u data, use i_shift_x = 0.0,  i_shift_y = -0.5
			 *             to interpolate given v data, use i_shift_y = -0.5, i_shift_y = 0.0
			 *  pay attention to the negative shift, which is necessary because the staggered grids are positively shifted
			 *  and this shift has to be removed for the interpolation
			 */
			double pos_x = wrapPeriodic(i_pos_x.scalar_data[pos_idx]*cached_scale_factor[0] + i_shift_x, (double)res[0]);
			double pos_y = wrapPeriodic(i_pos_y.scalar_data[pos_idx]*cached_scale_factor[1] + i_shift_y, (double)res[1]);

			/**
			 * For the interpolation, we assume node-aligned values
			 *
			 * x0  x1  x2  x3
			 * |---|---|---|---
			 * 0   2   4   6    <- positions and associated values e.g. for domain size 8
			 */

			/**
			 * See http://www.paulinternet.nl/?page=bicubic
			 */
			// compute x/y position
			double x = pos_x - floor(pos_x);
			double y = pos_y - floor(pos_y);

			// precompute x-position indices since they are reused 4 times
			int idx_i[10];  // up to 8th order center approximation to first derivative
			{

// TODO: Each 2nd wrapPeriodic is obsolete!!!
				int i = (int)pos_x-4; // C2
				/////////////int i = (int)pos_x-2;

				i = wrapPeriodic(i, res[0]);
				idx_i[0] = wrapPeriodic(i, res[0]);

                                for (int j = 1; j < 10; j++)
				{
					i = wrapPeriodic(i + 1, res[0]);
					idx_i[j] = wrapPeriodic(i, res[0]);
				}

			}

			/**
			 * iterate over rows and interpolate over the columns in the x direction
			 */
			// start at this row
			/////////////////int idx_j = wrapPeriodic((int)pos_y-2, res[1]);
			int idx_j = wrapPeriodic((int)pos_y-4, res[1]);

			double q[10];
			for (int kj = 0; kj < 10; kj++)
			{
				double p[10];

                                for (int j = 0; j < 10; j++)
					p[j] = i_data.physical_get(idx_j, idx_i[j]);

				double x0 = 0.;
				double x1 = 1.;
				double p0 = p[4];
				double p1 = p[5];

				double derx0 = 1. / 12. * p[2]  - 2. / 3. * p[3] + 2. / 3 * p[5] - 1. / 12 * p[6];
				double derx1 = 1. / 12. * p[3]  - 2. / 3. * p[4] + 2. / 3 * p[6] - 1. / 12 * p[7];

				double xx0 = x - x0;
				double xx0_2 = xx0 * xx0;
				double xx1 = x - x1;
				double x1x0 = x1 - x0;
	
				q[kj] = p0 + derx0 * xx0
					   + (p1 - p0 - derx0 * x1x0) * xx0_2 / (x1x0 * x1x0)
					   + ( (derx0 + derx1) * x1x0 - 2. * (p1 - p0) ) * xx0_2 * xx1 / (x1x0 * x1x0 * x1x0);

				idx_j = wrapPeriodic(idx_j+1, res[1]);
			}

			double y0 = 0.;
			double y1 = 1.;
			double q0 = q[4];
			double q1 = q[5];

			double dery0  = 1. / 12. * q[2]  - 2. / 3. * q[3] + 2. / 3 * q[5] - 1. / 12 * q[6];
			double dery1  = 1. / 12. * q[3]  - 2. / 3. * q[4] + 2. / 3 * q[6] - 1. / 12 * q[7];

			double yy0 = y - y0;
			double yy0_2 = yy0 * yy0;
			double yy1 = y - y1;
			double y1y0 = y1 - y0;

			double value = q0 + dery0 * yy0
					  + (q1 - q0 - dery0 * y1y0) * yy0_2 / (y1y0 * y1y0)
					  + ( (dery0 + dery1) * y1y0 - 2. * (q1 - q0) ) * yy0_2 * yy1 / (y1y0 * y1y0 * y1y0);

			o_data[pos_idx] = value;
		}
	}


public:
	void bicubic_scalar_NEW6(
			const PlaneData_Physical &i_data,			///< sampling data

			const ScalarDataArray &i_pos_x,		///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,		///< y positions of interpolation points

			double *o_data,						///< output values

			double i_shift_x = 0.0,				///< shift in x for staggered grids
			double i_shift_y = 0.0				///< shift in y for staggered grids
	)
	{
		assert(res[0] > 0);
		assert(cached_scale_factor[0] > 0);
		assert(i_pos_x.number_of_elements == i_pos_y.number_of_elements);

		std::size_t max_pos_idx = i_pos_x.number_of_elements;

#if SWEET_DEBUG
#if SWEET_THREADING_SPACE || SWEET_THREADING_TIME_REXI

		if (omp_get_num_threads() > 1)
			SWEETError("Are we already in parallel region? Threading race conditions likely!");
#endif
#endif

		// iterate over all positions in parallel
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t pos_idx = 0; pos_idx < max_pos_idx; pos_idx++)
		{
			/*
			 * load position to interpolate
			 * posx stores all x-coordinates of the arrival points
			 * posy stores all y-coordinates of the arrival points
			 *
			 * Both are arrays and matching array indices (pos_idx) below index the coordinates for the same point.
			 
			 * Scale factor (Nx/dx, Ny/dy) maps from the physical space to the array space.
			 * The array space is from [0; N[
			 *
			 * shift_x/y is operating in array space. Hence, staggered grid can be
			 * implemented by setting this to 0.5 or -0.5
			 * for C grid, to interpolate given u data, use i_shift_x = 0.0,  i_shift_y = -0.5
			 *             to interpolate given v data, use i_shift_y = -0.5, i_shift_y = 0.0
			 *  pay attention to the negative shift, which is necessary because the staggered grids are positively shifted
			 *  and this shift has to be removed for the interpolation
			 */
			double pos_x = wrapPeriodic(i_pos_x.scalar_data[pos_idx]*cached_scale_factor[0] + i_shift_x, (double)res[0]);
			double pos_y = wrapPeriodic(i_pos_y.scalar_data[pos_idx]*cached_scale_factor[1] + i_shift_y, (double)res[1]);

			/**
			 * For the interpolation, we assume node-aligned values
			 *
			 * x0  x1  x2  x3
			 * |---|---|---|---
			 * 0   2   4   6    <- positions and associated values e.g. for domain size 8
			 */

			/**
			 * See http://www.paulinternet.nl/?page=bicubic
			 */
			// compute x/y position
			double x = pos_x - floor(pos_x);
			double y = pos_y - floor(pos_y);

			// precompute x-position indices since they are reused 4 times
			int idx_i[10];  // up to 8th order center approximation to first derivative
			{

// TODO: Each 2nd wrapPeriodic is obsolete!!!
				int i = (int)pos_x-4; // C2
				/////////////int i = (int)pos_x-2;

				i = wrapPeriodic(i, res[0]);
				idx_i[0] = wrapPeriodic(i, res[0]);

                                for (int j = 1; j < 10; j++)
				{
					i = wrapPeriodic(i + 1, res[0]);
					idx_i[j] = wrapPeriodic(i, res[0]);
				}

			}

			/**
			 * iterate over rows and interpolate over the columns in the x direction
			 */
			// start at this row
			/////////////////int idx_j = wrapPeriodic((int)pos_y-2, res[1]);
			int idx_j = wrapPeriodic((int)pos_y-4, res[1]);

			double q[10];
			for (int kj = 0; kj < 10; kj++)
			{
				double p[10];

                                for (int j = 0; j < 10; j++)
					p[j] = i_data.physical_get(idx_j, idx_i[j]);

				double x0 = 0.;
				double x1 = 1.;
				double p0 = p[4];
				double p1 = p[5];

				double derx0 = -1. / 60. * p[1] + 3. / 20. * p[2] - 3. / 4. * p[3] + 3. / 4. * p[5] - 3. / 20. * p[6] + 1. / 60. * p[7];
				double derx1 = -1. / 60. * p[2] + 3. / 20. * p[3] - 3. / 4. * p[4] + 3. / 4. * p[6] - 3. / 20. * p[7] + 1. / 60. * p[8];

				double xx0 = x - x0;
				double xx0_2 = xx0 * xx0;
				double xx1 = x - x1;
				double x1x0 = x1 - x0;
	
				q[kj] = p0 + derx0 * xx0
					   + (p1 - p0 - derx0 * x1x0) * xx0_2 / (x1x0 * x1x0)
					   + ( (derx0 + derx1) * x1x0 - 2. * (p1 - p0) ) * xx0_2 * xx1 / (x1x0 * x1x0 * x1x0);

				idx_j = wrapPeriodic(idx_j+1, res[1]);
			}

			double y0 = 0.;
			double y1 = 1.;
			double q0 = q[4];
			double q1 = q[5];

			double dery0  = -1. / 60. * q[1] + 3. / 20. * q[2] - 3. / 4. * q[3] + 3. / 4. * q[5] - 3. / 20. * q[6] + 1. / 60. * q[7];
			double dery1  = -1. / 60. * q[2] + 3. / 20. * q[3] - 3. / 4. * q[4] + 3. / 4. * q[6] - 3. / 20. * q[7] + 1. / 60. * q[8];

			double yy0 = y - y0;
			double yy0_2 = yy0 * yy0;
			double yy1 = y - y1;
			double y1y0 = y1 - y0;

			double value = q0 + dery0 * yy0
					  + (q1 - q0 - dery0 * y1y0) * yy0_2 / (y1y0 * y1y0)
					  + ( (dery0 + dery1) * y1y0 - 2. * (q1 - q0) ) * yy0_2 * yy1 / (y1y0 * y1y0 * y1y0);

			o_data[pos_idx] = value;
		}
	}


public:
	void bicubic_scalar_NEW8(
			const PlaneData_Physical &i_data,			///< sampling data

			const ScalarDataArray &i_pos_x,		///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,		///< y positions of interpolation points

			double *o_data,						///< output values

			double i_shift_x = 0.0,				///< shift in x for staggered grids
			double i_shift_y = 0.0				///< shift in y for staggered grids
	)
	{
		assert(res[0] > 0);
		assert(cached_scale_factor[0] > 0);
		assert(i_pos_x.number_of_elements == i_pos_y.number_of_elements);

		std::size_t max_pos_idx = i_pos_x.number_of_elements;

#if SWEET_DEBUG
#if SWEET_THREADING_SPACE || SWEET_THREADING_TIME_REXI

		if (omp_get_num_threads() > 1)
			SWEETError("Are we already in parallel region? Threading race conditions likely!");
#endif
#endif

		// iterate over all positions in parallel
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t pos_idx = 0; pos_idx < max_pos_idx; pos_idx++)
		{
			/*
			 * load position to interpolate
			 * posx stores all x-coordinates of the arrival points
			 * posy stores all y-coordinates of the arrival points
			 *
			 * Both are arrays and matching array indices (pos_idx) below index the coordinates for the same point.
			 
			 * Scale factor (Nx/dx, Ny/dy) maps from the physical space to the array space.
			 * The array space is from [0; N[
			 *
			 * shift_x/y is operating in array space. Hence, staggered grid can be
			 * implemented by setting this to 0.5 or -0.5
			 * for C grid, to interpolate given u data, use i_shift_x = 0.0,  i_shift_y = -0.5
			 *             to interpolate given v data, use i_shift_y = -0.5, i_shift_y = 0.0
			 *  pay attention to the negative shift, which is necessary because the staggered grids are positively shifted
			 *  and this shift has to be removed for the interpolation
			 */
			double pos_x = wrapPeriodic(i_pos_x.scalar_data[pos_idx]*cached_scale_factor[0] + i_shift_x, (double)res[0]);
			double pos_y = wrapPeriodic(i_pos_y.scalar_data[pos_idx]*cached_scale_factor[1] + i_shift_y, (double)res[1]);

			/**
			 * For the interpolation, we assume node-aligned values
			 *
			 * x0  x1  x2  x3
			 * |---|---|---|---
			 * 0   2   4   6    <- positions and associated values e.g. for domain size 8
			 */

			/**
			 * See http://www.paulinternet.nl/?page=bicubic
			 */
			// compute x/y position
			double x = pos_x - floor(pos_x);
			double y = pos_y - floor(pos_y);

			// precompute x-position indices since they are reused 4 times
			int idx_i[10];  // up to 8th order center approximation to first derivative
			{

// TODO: Each 2nd wrapPeriodic is obsolete!!!
				int i = (int)pos_x-4; // C2
				/////////////int i = (int)pos_x-2;

				i = wrapPeriodic(i, res[0]);
				idx_i[0] = wrapPeriodic(i, res[0]);

                                for (int j = 1; j < 10; j++)
				{
					i = wrapPeriodic(i + 1, res[0]);
					idx_i[j] = wrapPeriodic(i, res[0]);
				}

			}

			/**
			 * iterate over rows and interpolate over the columns in the x direction
			 */
			// start at this row
			/////////////////int idx_j = wrapPeriodic((int)pos_y-2, res[1]);
			int idx_j = wrapPeriodic((int)pos_y-4, res[1]);

			double q[10];
			for (int kj = 0; kj < 10; kj++)
			{
				double p[10];

                                for (int j = 0; j < 10; j++)
					p[j] = i_data.physical_get(idx_j, idx_i[j]);

				double x0 = 0.;
				double x1 = 1.;
				double p0 = p[4];
				double p1 = p[5];

				double derx0  = 1. / 280. * p[0] - 4. / 105. * p[1] + 1. / 5. * p[2] - 4. / 5 * p[3] + 4. / 5 * p[5] - 1. / 5. * p[6] + 4. / 105. * p[7] - 1. / 280. * p[8];
				double derx1  = 1. / 280. * p[1] - 4. / 105. * p[2] + 1. / 5. * p[3] - 4. / 5 * p[4] + 4. / 5 * p[6] - 1. / 5. * p[7] + 4. / 105. * p[8] - 1. / 280. * p[9];

				double xx0 = x - x0;
				double xx0_2 = xx0 * xx0;
				double xx1 = x - x1;
				double x1x0 = x1 - x0;
	
				q[kj] = p0 + derx0 * xx0
					   + (p1 - p0 - derx0 * x1x0) * xx0_2 / (x1x0 * x1x0)
					   + ( (derx0 + derx1) * x1x0 - 2. * (p1 - p0) ) * xx0_2 * xx1 / (x1x0 * x1x0 * x1x0);

				idx_j = wrapPeriodic(idx_j+1, res[1]);
			}

			double y0 = 0.;
			double y1 = 1.;
			double q0 = q[4];
			double q1 = q[5];

			double dery0  = 1. / 280. * q[0] - 4. / 105. * q[1] + 1. / 5. * q[2] - 4. / 5 * q[3] + 4. / 5 * q[5] - 1. / 5. * q[6] + 4. / 105. * q[7] - 1. / 280. * q[8];
			double dery1  = 1. / 280. * q[1] - 4. / 105. * q[2] + 1. / 5. * q[3] - 4. / 5 * q[4] + 4. / 5 * q[6] - 1. / 5. * q[7] + 4. / 105. * q[8] - 1. / 280. * q[9];

			double yy0 = y - y0;
			double yy0_2 = yy0 * yy0;
			double yy1 = y - y1;
			double y1y0 = y1 - y0;

			double value = q0 + dery0 * yy0
					  + (q1 - q0 - dery0 * y1y0) * yy0_2 / (y1y0 * y1y0)
					  + ( (dery0 + dery1) * y1y0 - 2. * (q1 - q0) ) * yy0_2 * yy1 / (y1y0 * y1y0 * y1y0);

			o_data[pos_idx] = value;
		}
	}

public:
	void bicubic_scalar_Lagrangian(
			const PlaneData_Physical &i_data,			///< sampling data

			const ScalarDataArray &i_pos_x,		///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,		///< y positions of interpolation points

			double *o_data,						///< output values

			double i_shift_x = 0.0,				///< shift in x for staggered grids
			double i_shift_y = 0.0				///< shift in y for staggered grids
	)
	{
		assert(res[0] > 0);
		assert(cached_scale_factor[0] > 0);
		assert(i_pos_x.number_of_elements == i_pos_y.number_of_elements);

		std::size_t max_pos_idx = i_pos_x.number_of_elements;

#if SWEET_DEBUG
#if SWEET_THREADING_SPACE || SWEET_THREADING_TIME_REXI

		if (omp_get_num_threads() > 1)
			SWEETError("Are we already in parallel region? Threading race conditions likely!");
#endif
#endif

		// iterate over all positions in parallel
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t pos_idx = 0; pos_idx < max_pos_idx; pos_idx++)
		{
			/*
			 * load position to interpolate
			 * posx stores all x-coordinates of the arrival points
			 * posy stores all y-coordinates of the arrival points
			 *
			 * Both are arrays and matching array indices (pos_idx) below index the coordinates for the same point.
			 
			 * Scale factor (Nx/dx, Ny/dy) maps from the physical space to the array space.
			 * The array space is from [0; N[
			 *
			 * shift_x/y is operating in array space. Hence, staggered grid can be
			 * implemented by setting this to 0.5 or -0.5
			 * for C grid, to interpolate given u data, use i_shift_x = 0.0,  i_shift_y = -0.5
			 *             to interpolate given v data, use i_shift_y = -0.5, i_shift_y = 0.0
			 *  pay attention to the negative shift, which is necessary because the staggered grids are positively shifted
			 *  and this shift has to be removed for the interpolation
			 */
			double pos_x = wrapPeriodic(i_pos_x.scalar_data[pos_idx]*cached_scale_factor[0] + i_shift_x, (double)res[0]);
			double pos_y = wrapPeriodic(i_pos_y.scalar_data[pos_idx]*cached_scale_factor[1] + i_shift_y, (double)res[1]);

			/**
			 * For the interpolation, we assume node-aligned values
			 *
			 * x0  x1  x2  x3
			 * |---|---|---|---
			 * 0   2   4   6    <- positions and associated values e.g. for domain size 8
			 */

			/**
			 * See http://www.paulinternet.nl/?page=bicubic
			 */
			// compute x/y position
			double x = pos_x - floor(pos_x);
			double y = pos_y - floor(pos_y);

			// precompute x-position indices since they are reused 4 times
			int idx_i[10];  // up to 8th order center approximation to first derivative
			{

// TODO: Each 2nd wrapPeriodic is obsolete!!!
				int i = (int)pos_x-4; // C2
				/////////////int i = (int)pos_x-2;

				i = wrapPeriodic(i, res[0]);
				idx_i[0] = wrapPeriodic(i, res[0]);

                                for (int j = 1; j < 10; j++)
				{
					i = wrapPeriodic(i + 1, res[0]);
					idx_i[j] = wrapPeriodic(i, res[0]);
				}

			}

			/**
			 * iterate over rows and interpolate over the columns in the x direction
			 */
			// start at this row
			/////////////////int idx_j = wrapPeriodic((int)pos_y-2, res[1]);
			int idx_j = wrapPeriodic((int)pos_y-4, res[1]);

			double q[10];
			for (int kj = 0; kj < 10; kj++)
			{
				double p[10];

                                for (int j = 0; j < 10; j++)
					p[j] = i_data.physical_get(idx_j, idx_i[j]);

				double x0 = -1.;
				double x1 = 0.;
				double x2 = 1.;
				double x3 = 2.;
				double p0 = p[3];
				double p1 = p[4];
				double p2 = p[5];
				double p3 = p[6];

				double xx0 = x - x0;
				double xx1 = (x - x1) * xx0;
				double xx2 = (x - x2) * xx1;
	
				q[kj] = p0 + (p1 - p0) * xx0
					   + (p2 - 2. * p1 + p0) / 2. * xx1
					   + (p3 - 3. * p2 + 3. * p1 - p0) / 6. * xx2;


				idx_j = wrapPeriodic(idx_j+1, res[1]);
			}

			double y0 = -1.;
			double y1 = 0.;
			double y2 = 1.;
			double y3 = 2.;
			double q0 = q[3];
			double q1 = q[4];
			double q2 = q[5];
			double q3 = q[6];

			double yy0 = y - y0;
			double yy1 = (y - y1) * yy0;
			double yy2 = (y - y2) * yy1;

			double value = q0 + (q1 - q0) * yy0
                                          + (q2 - 2. * q1 + q0) / 2. * yy1
                                          + (q3 - 3. * q2 + 3. * q1 - q0) / 6. * yy2;

			o_data[pos_idx] = value;
		}
	}



public:
	void biquintic_scalar_Lagrangian(
			const PlaneData_Physical &i_data,			///< sampling data

			const ScalarDataArray &i_pos_x,		///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,		///< y positions of interpolation points

			double *o_data,						///< output values

			double i_shift_x = 0.0,				///< shift in x for staggered grids
			double i_shift_y = 0.0				///< shift in y for staggered grids
	)
	{
		assert(res[0] > 0);
		assert(cached_scale_factor[0] > 0);
		assert(i_pos_x.number_of_elements == i_pos_y.number_of_elements);

		std::size_t max_pos_idx = i_pos_x.number_of_elements;

#if SWEET_DEBUG
#if SWEET_THREADING_SPACE || SWEET_THREADING_TIME_REXI

		if (omp_get_num_threads() > 1)
			SWEETError("Are we already in parallel region? Threading race conditions likely!");
#endif
#endif

		// iterate over all positions in parallel
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t pos_idx = 0; pos_idx < max_pos_idx; pos_idx++)
		{
			/*
			 * load position to interpolate
			 * posx stores all x-coordinates of the arrival points
			 * posy stores all y-coordinates of the arrival points
			 *
			 * Both are arrays and matching array indices (pos_idx) below index the coordinates for the same point.
			 
			 * Scale factor (Nx/dx, Ny/dy) maps from the physical space to the array space.
			 * The array space is from [0; N[
			 *
			 * shift_x/y is operating in array space. Hence, staggered grid can be
			 * implemented by setting this to 0.5 or -0.5
			 * for C grid, to interpolate given u data, use i_shift_x = 0.0,  i_shift_y = -0.5
			 *             to interpolate given v data, use i_shift_y = -0.5, i_shift_y = 0.0
			 *  pay attention to the negative shift, which is necessary because the staggered grids are positively shifted
			 *  and this shift has to be removed for the interpolation
			 */
			double pos_x = wrapPeriodic(i_pos_x.scalar_data[pos_idx]*cached_scale_factor[0] + i_shift_x, (double)res[0]);
			double pos_y = wrapPeriodic(i_pos_y.scalar_data[pos_idx]*cached_scale_factor[1] + i_shift_y, (double)res[1]);

			/**
			 * For the interpolation, we assume node-aligned values
			 *
			 * x0  x1  x2  x3
			 * |---|---|---|---
			 * 0   2   4   6    <- positions and associated values e.g. for domain size 8
			 */

			/**
			 * See http://www.paulinternet.nl/?page=bicubic
			 */
			// compute x/y position
			double x = pos_x - floor(pos_x);
			double y = pos_y - floor(pos_y);

			// precompute x-position indices since they are reused 4 times
			int idx_i[10];  // up to 8th order center approximation to first derivative
			{

// TODO: Each 2nd wrapPeriodic is obsolete!!!
				int i = (int)pos_x-4; // C2
				/////////////int i = (int)pos_x-2;

				i = wrapPeriodic(i, res[0]);
				idx_i[0] = wrapPeriodic(i, res[0]);

                                for (int j = 1; j < 10; j++)
				{
					i = wrapPeriodic(i + 1, res[0]);
					idx_i[j] = wrapPeriodic(i, res[0]);
				}

			}

			/**
			 * iterate over rows and interpolate over the columns in the x direction
			 */
			// start at this row
			/////////////////int idx_j = wrapPeriodic((int)pos_y-2, res[1]);
			int idx_j = wrapPeriodic((int)pos_y-4, res[1]);

			double q[10];
			for (int kj = 0; kj < 10; kj++)
			{
				double p[10];

                                for (int j = 0; j < 10; j++)
					p[j] = i_data.physical_get(idx_j, idx_i[j]);

				double x0 = -2.;
				double x1 = -1.;
				double x2 = 0.;
				double x3 = 1.;
				double x4 = 2.;
				double x5 = 3.;
				double p0 = p[2];
				double p1 = p[3];
				double p2 = p[4];
				double p3 = p[5];
				double p4 = p[6];
				double p5 = p[7];

				double xx0 = x - x0;
				double xx1 = (x - x1) * xx0;
				double xx2 = (x - x2) * xx1;
				double xx3 = (x - x3) * xx2;
				double xx4 = (x - x4) * xx3;
	
				q[kj] = p0 + (p1 - p0) * xx0
					   + (p2 - 2. * p1 + p0) / 2. * xx1
					   + (p3 - 3. * p2 + 3. * p1 - p0) / 6. * xx2
					   + (p4 - 4. * p3 + 6. * p2 - 4. * p1 + p0) / 24. * xx3
					   + (p5 - 5. * p4 + 10. * p3 - 10. * p2 + 5. * p1 - p0) / 120. * xx4;


				idx_j = wrapPeriodic(idx_j+1, res[1]);
			}

			double y0 = -2.;
			double y1 = -1.;
			double y2 = 0.;
			double y3 = 1.;
			double y4 = 2.;
			double y5 = 3.;
			double q0 = q[2];
			double q1 = q[3];
			double q2 = q[4];
			double q3 = q[5];
			double q4 = q[6];
			double q5 = q[7];

			double yy0 = y - y0;
			double yy1 = (y - y1) * yy0;
			double yy2 = (y - y2) * yy1;
			double yy3 = (y - y3) * yy2;
			double yy4 = (y - y4) * yy3;

			double value = q0 + (q1 - q0) * yy0
                                          + (q2 - 2. * q1 + q0) / 2. * yy1
                                          + (q3 - 3. * q2 + 3. * q1 - q0) / 6. * yy2
					  + (q4 - 4. * q3 + 6. * q2 - 4. * q1 + q0) / 24. * yy3
					  + (q5 - 5. * q4 + 10. * q3 - 10. * q2 + 5. * q1 - q0) / 120. * yy4;

			o_data[pos_idx] = value;
		}
	}


public:
	void biseptic_scalar_Lagrangian(
			const PlaneData_Physical &i_data,			///< sampling data

			const ScalarDataArray &i_pos_x,		///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,		///< y positions of interpolation points

			double *o_data,						///< output values

			double i_shift_x = 0.0,				///< shift in x for staggered grids
			double i_shift_y = 0.0				///< shift in y for staggered grids
	)
	{
		assert(res[0] > 0);
		assert(cached_scale_factor[0] > 0);
		assert(i_pos_x.number_of_elements == i_pos_y.number_of_elements);

		std::size_t max_pos_idx = i_pos_x.number_of_elements;

#if SWEET_DEBUG
#if SWEET_THREADING_SPACE || SWEET_THREADING_TIME_REXI

		if (omp_get_num_threads() > 1)
			SWEETError("Are we already in parallel region? Threading race conditions likely!");
#endif
#endif

		// iterate over all positions in parallel
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t pos_idx = 0; pos_idx < max_pos_idx; pos_idx++)
		{
			/*
			 * load position to interpolate
			 * posx stores all x-coordinates of the arrival points
			 * posy stores all y-coordinates of the arrival points
			 *
			 * Both are arrays and matching array indices (pos_idx) below index the coordinates for the same point.
			 
			 * Scale factor (Nx/dx, Ny/dy) maps from the physical space to the array space.
			 * The array space is from [0; N[
			 *
			 * shift_x/y is operating in array space. Hence, staggered grid can be
			 * implemented by setting this to 0.5 or -0.5
			 * for C grid, to interpolate given u data, use i_shift_x = 0.0,  i_shift_y = -0.5
			 *             to interpolate given v data, use i_shift_y = -0.5, i_shift_y = 0.0
			 *  pay attention to the negative shift, which is necessary because the staggered grids are positively shifted
			 *  and this shift has to be removed for the interpolation
			 */
			double pos_x = wrapPeriodic(i_pos_x.scalar_data[pos_idx]*cached_scale_factor[0] + i_shift_x, (double)res[0]);
			double pos_y = wrapPeriodic(i_pos_y.scalar_data[pos_idx]*cached_scale_factor[1] + i_shift_y, (double)res[1]);

			/**
			 * For the interpolation, we assume node-aligned values
			 *
			 * x0  x1  x2  x3
			 * |---|---|---|---
			 * 0   2   4   6    <- positions and associated values e.g. for domain size 8
			 */

			/**
			 * See http://www.paulinternet.nl/?page=bicubic
			 */
			// compute x/y position
			double x = pos_x - floor(pos_x);
			double y = pos_y - floor(pos_y);

			// precompute x-position indices since they are reused 4 times
			int idx_i[10];  // up to 8th order center approximation to first derivative
			{

// TODO: Each 2nd wrapPeriodic is obsolete!!!
				int i = (int)pos_x-4; // C2
				/////////////int i = (int)pos_x-2;

				i = wrapPeriodic(i, res[0]);
				idx_i[0] = wrapPeriodic(i, res[0]);

                                for (int j = 1; j < 10; j++)
				{
					i = wrapPeriodic(i + 1, res[0]);
					idx_i[j] = wrapPeriodic(i, res[0]);
				}

			}

			/**
			 * iterate over rows and interpolate over the columns in the x direction
			 */
			// start at this row
			/////////////////int idx_j = wrapPeriodic((int)pos_y-2, res[1]);
			int idx_j = wrapPeriodic((int)pos_y-4, res[1]);

			double q[10];
			for (int kj = 0; kj < 10; kj++)
			{
				double p[10];

                                for (int j = 0; j < 10; j++)
					p[j] = i_data.physical_get(idx_j, idx_i[j]);

				double x0 = -3.;
				double x1 = -2.;
				double x2 = -1.;
				double x3 = 0.;
				double x4 = 1.;
				double x5 = 2.;
				double x6 = 3.;
				double x7 = 4.;
				double p0 = p[1];
				double p1 = p[2];
				double p2 = p[3];
				double p3 = p[4];
				double p4 = p[5];
				double p5 = p[6];
				double p6 = p[7];
				double p7 = p[8];

				double xx0 = x - x0;
				double xx1 = (x - x1) * xx0;
				double xx2 = (x - x2) * xx1;
				double xx3 = (x - x3) * xx2;
				double xx4 = (x - x4) * xx3;
				double xx5 = (x - x5) * xx4;
				double xx6 = (x - x6) * xx5;
	
				q[kj] = p0 + (p1 - p0) * xx0
					   + (p2 - 2. * p1 + p0) / 2. * xx1
					   + (p3 - 3. * p2 + 3. * p1 - p0) / 6. * xx2
					   + (p4 - 4. * p3 + 6. * p2 - 4. * p1 + p0) / 24. * xx3
					   + (p5 - 5. * p4 + 10. * p3 - 10. * p2 + 5. * p1 - p0) / 120. * xx4
					   + (p6 - 6. * p5 + 15. * p4 - 20. * p3 + 15. * p2 - 6. * p1 + p0) / 720. * xx5
					   + (p7 - 7. * p6 + 21. * p5 - 35. * p4 + 35. * p3 - 21. * p2 + 7. * p1 - p0) / 5040. * xx6;


				idx_j = wrapPeriodic(idx_j+1, res[1]);
			}

			double y0 = -3.;
			double y1 = -2.;
			double y2 = -1.;
			double y3 = 0.;
			double y4 = 1.;
			double y5 = 2.;
			double y6 = 3.;
			double y7 = 4.;
			double q0 = q[1];
			double q1 = q[2];
			double q2 = q[3];
			double q3 = q[4];
			double q4 = q[5];
			double q5 = q[6];
			double q6 = q[7];
			double q7 = q[8];

			double yy0 = y - y0;
			double yy1 = (y - y1) * yy0;
			double yy2 = (y - y2) * yy1;
			double yy3 = (y - y3) * yy2;
			double yy4 = (y - y4) * yy3;
			double yy5 = (y - y5) * yy4;
			double yy6 = (y - y6) * yy5;

			double value = q0 + (q1 - q0) * yy0
                                          + (q2 - 2. * q1 + q0) / 2. * yy1
                                          + (q3 - 3. * q2 + 3. * q1 - q0) / 6. * yy2
					  + (q4 - 4. * q3 + 6. * q2 - 4. * q1 + q0) / 24. * yy3
					  + (q5 - 5. * q4 + 10. * q3 - 10. * q2 + 5. * q1 - q0) / 120. * yy4
					  + (q6 - 6. * q5 + 15. * q4 - 20. * q3 + 15. * q2 - 6. * q1 + q0) / 720. * yy5
					  + (q7 - 7. * q6 + 21. * q5 - 35. * q4 + 35. * q3 - 21. * q2 + 7. * q1 - q0) / 5040. * yy6;

			o_data[pos_idx] = value;
		}
	}


public:
	void biquintic_scalar(
			const PlaneData_Physical &i_data,			///< sampling data

			const ScalarDataArray &i_pos_x,		///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,		///< y positions of interpolation points

			double *o_data,						///< output values

			double i_shift_x = 0.0,				///< shift in x for staggered grids
			double i_shift_y = 0.0				///< shift in y for staggered grids
	)
	{
		assert(res[0] > 0);
		assert(cached_scale_factor[0] > 0);
		assert(i_pos_x.number_of_elements == i_pos_y.number_of_elements);

		std::size_t max_pos_idx = i_pos_x.number_of_elements;

#if SWEET_DEBUG
#if SWEET_THREADING_SPACE || SWEET_THREADING_TIME_REXI

		if (omp_get_num_threads() > 1)
			SWEETError("Are we already in parallel region? Threading race conditions likely!");
#endif
#endif

		// iterate over all positions in parallel
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t pos_idx = 0; pos_idx < max_pos_idx; pos_idx++)
		{
			/*
			 * load position to interpolate
			 * posx stores all x-coordinates of the arrival points
			 * posy stores all y-coordinates of the arrival points
			 *
			 * Both are arrays and matching array indices (pos_idx) below index the coordinates for the same point.
			 *
			 * Scale factor (Nx/dx, Ny/dy) maps from the physical space to the array space.
			 * The array space is from [0; N[
			 *
			 * shift_x/y is operating in array space. Hence, staggered grid can be
			 * implemented by setting this to 0.5 or -0.5
			 * for C grid, to interpolate given u data, use i_shift_x = 0.0,  i_shift_y = -0.5
			 *             to interpolate given v data, use i_shift_y = -0.5, i_shift_y = 0.0
			 *  pay attention to the negative shift, which is necessary because the staggered grids are positively shifted
			 *  and this shift has to be removed for the interpolation
			 */
			double pos_x = wrapPeriodic(i_pos_x.scalar_data[pos_idx]*cached_scale_factor[0] + i_shift_x, (double)res[0]);
			double pos_y = wrapPeriodic(i_pos_y.scalar_data[pos_idx]*cached_scale_factor[1] + i_shift_y, (double)res[1]);

			/**
			 * For the interpolation, we assume node-aligned values
			 *
			 * x0  x1  x2  x3
			 * |---|---|---|---
			 * 0   2   4   6    <- positions and associated values e.g. for domain size 8
			 */

			/**
			 * See http://www.paulinternet.nl/?page=bicubic
			 */
			// compute x/y position
			double x = pos_x - floor(pos_x);
			double y = pos_y - floor(pos_y);

			// precompute x-position indices since they are reused 4 times
			int idx_i[6];
			{

// TODO: Each 2nd wrapPeriodic is obsolete!!!
				int i = (int)pos_x-1;
				/////////////int i = (int)pos_x-2;

				i = wrapPeriodic(i, res[0]);
				idx_i[0] = wrapPeriodic(i, res[0]);

				i = wrapPeriodic(i+1, res[0]);
				idx_i[1] = wrapPeriodic(i, res[0]);

				i = wrapPeriodic(i+1, res[0]);
				idx_i[2] = wrapPeriodic(i, res[0]);

				i = wrapPeriodic(i+1, res[0]);
				idx_i[3] = wrapPeriodic(i, res[0]);

				i = wrapPeriodic(i+1, res[0]);
				idx_i[4] = wrapPeriodic(i, res[0]);

				i = wrapPeriodic(i+1, res[0]);
				idx_i[5] = wrapPeriodic(i, res[0]);
			}

			/**
			 * iterate over rows and interpolate over the columns in the x direction
			 */
			// start at this row
			/////////////////int idx_j = wrapPeriodic((int)pos_y-2, res[1]);
			int idx_j = wrapPeriodic((int)pos_y-1, res[1]);

			double q[6];
			for (int kj = 0; kj < 6; kj++)
			{
				double p[6];

				p[0] = i_data.physical_get(idx_j, idx_i[0]);
				p[1] = i_data.physical_get(idx_j, idx_i[1]);
				p[2] = i_data.physical_get(idx_j, idx_i[2]);
				p[3] = i_data.physical_get(idx_j, idx_i[3]);
				p[4] = i_data.physical_get(idx_j, idx_i[4]);
				p[5] = i_data.physical_get(idx_j, idx_i[5]);

				double x0 = 0.;
				double x1 = 1.;
				double p0 = p[1];
				double p1 = p[2];

				double derx0  = (p[2] - p[0]) / 2.;
				double derx1  = (p[3] - p[1]) / 2.;
				double derxx0 = p[2] - 2. * p[1] + p[0];
				double derxx1 = p[3] - 2. * p[2] + p[1];

				double xx0 = x - x0;
				double xx0_2 = xx0 * xx0;
				double xx0_3 = xx0_2 * xx0;
				double xx1 = x - x1;
				double xx1_2 = xx1 * xx1;
				double xx1_3 = xx1_2 * xx1;
				double x1x0 = x1 -x0;
				double x1x0_2 = x1x0 * x1x0;
				double x1x0_3 = x1x0_2 * x1x0;
				double x1x0_4 = x1x0_3 * x1x0;
				double x1x0_5 = x1x0_4 * x1x0;
	
				q[kj] = p0 + derx0 * xx0 + .5 * derxx0 * xx0_2
				        + ( p1 - p0 - derx0 * x1x0 - .5 * derxx0 * x1x0_2) * xx0_3 / x1x0_3
				        + ( 3. * p0 - 3. * p1  + (2. * derx0 + derx1) * x1x0 + .5 * derxx0 * x1x0_2  ) * xx0_3 * xx1 / x1x0_4
				        + ( 6. * p1 - 6. * p0 - 3. * (derx0 + derx1) * x1x0 + .5 * (derxx1 - derxx0) * x1x0_2  ) * xx0_3 * xx1_2 / x1x0_5;

				idx_j = wrapPeriodic(idx_j+1, res[1]);
			}

			double y0 = 0.;
			double y1 = 1.;
			double q0 = q[1];
			double q1 = q[2];

			double dery0  = (q[2] - q[0]) / 2.;
			double dery1  = (q[3] - q[1]) / 2.;
			double deryy0 = q[2] - 2. * q[1] + q[0];
			double deryy1 = q[3] - 2. * q[2] + q[1];

			double yy0 = y - y0;
			double yy0_2 = yy0 * yy0;
			double yy0_3 = yy0_2 * yy0;
			double yy1 = y - y1;
			double yy1_2 = yy1 * yy1;
			double yy1_3 = yy1_2 * yy1;
			double y1y0 = y1 -y0;
			double y1y0_2 = y1y0 * y1y0;
			double y1y0_3 = y1y0_2 * y1y0;
			double y1y0_4 = y1y0_3 * y1y0;
			double y1y0_5 = y1y0_4 * y1y0;

			double value = q0 + dery0 * yy0 + .5 * deryy0 * yy0_2
			               + ( q1 - q0 - dery0 * y1y0 - .5 * deryy0 * y1y0_2) * yy0_3 / y1y0_3
			               + ( 3. * q0 - 3. * q1  + (2. * dery0 + dery1) * y1y0 + .5 * deryy0 * y1y0_2  ) * yy0_3 * yy1 / y1y0_4
			               + ( 6. * q1 - 6. * q0 - 3. * (dery0 + dery1) * y1y0 + .5 * (deryy1 - deryy0) * y1y0_2  ) * yy0_3 * yy1_2 / y1y0_5;

			o_data[pos_idx] = value;
		}
	}




public:
	void biquintic_scalar2(
			const PlaneData_Physical &i_data,			///< sampling data

			const ScalarDataArray &i_pos_x,		///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,		///< y positions of interpolation points

			double *o_data,						///< output values

			double i_shift_x = 0.0,				///< shift in x for staggered grids
			double i_shift_y = 0.0				///< shift in y for staggered grids
	)
	{
		assert(res[0] > 0);
		assert(cached_scale_factor[0] > 0);
		assert(i_pos_x.number_of_elements == i_pos_y.number_of_elements);

		std::size_t max_pos_idx = i_pos_x.number_of_elements;

#if SWEET_DEBUG
#if SWEET_THREADING_SPACE || SWEET_THREADING_TIME_REXI

		if (omp_get_num_threads() > 1)
			SWEETError("Are we already in parallel region? Threading race conditions likely!");
#endif
#endif

		// iterate over all positions in parallel
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t pos_idx = 0; pos_idx < max_pos_idx; pos_idx++)
		{
			/*
			 * load position to interpolate
			 * posx stores all x-coordinates of the arrival points
			 * posy stores all y-coordinates of the arrival points
			 *
			 * Both are arrays and matching array indices (pos_idx) below index the coordinates for the same point.
			 *
			 * Scale factor (Nx/dx, Ny/dy) maps from the physical space to the array space.
			 * The array space is from [0; N[
			 *
			 * shift_x/y is operating in array space. Hence, staggered grid can be
			 * implemented by setting this to 0.5 or -0.5
			 * for C grid, to interpolate given u data, use i_shift_x = 0.0,  i_shift_y = -0.5
			 *             to interpolate given v data, use i_shift_y = -0.5, i_shift_y = 0.0
			 *  pay attention to the negative shift, which is necessary because the staggered grids are positively shifted
			 *  and this shift has to be removed for the interpolation
			 */
			double pos_x = wrapPeriodic(i_pos_x.scalar_data[pos_idx]*cached_scale_factor[0] + i_shift_x, (double)res[0]);
			double pos_y = wrapPeriodic(i_pos_y.scalar_data[pos_idx]*cached_scale_factor[1] + i_shift_y, (double)res[1]);

			/**
			 * For the interpolation, we assume node-aligned values
			 *
			 * x0  x1  x2  x3
			 * |---|---|---|---
			 * 0   2   4   6    <- positions and associated values e.g. for domain size 8
			 */

			/**
			 * See http://www.paulinternet.nl/?page=bicubic
			 */
			// compute x/y position
			double x = pos_x - floor(pos_x);
			double y = pos_y - floor(pos_y);

			// precompute x-position indices since they are reused 4 times
			int idx_i[6];
			{

// TODO: Each 2nd wrapPeriodic is obsolete!!!
				int i = (int)pos_x-2;

				i = wrapPeriodic(i, res[0]);
				idx_i[0] = wrapPeriodic(i, res[0]);

				i = wrapPeriodic(i+1, res[0]);
				idx_i[1] = wrapPeriodic(i, res[0]);

				i = wrapPeriodic(i+1, res[0]);
				idx_i[2] = wrapPeriodic(i, res[0]);

				i = wrapPeriodic(i+1, res[0]);
				idx_i[3] = wrapPeriodic(i, res[0]);

				i = wrapPeriodic(i+1, res[0]);
				idx_i[4] = wrapPeriodic(i, res[0]);

				i = wrapPeriodic(i+1, res[0]);
				idx_i[5] = wrapPeriodic(i, res[0]);
			}

			/**
			 * iterate over rows and interpolate over the columns in the x direction
			 */
			// start at this row
			/////////////////int idx_j = wrapPeriodic((int)pos_y-2, res[1]);
			int idx_j = wrapPeriodic((int)pos_y - 2, res[1]);

			double q[6];
			for (int kj = 0; kj < 6; kj++)
			{
				double p[6];

				p[0] = i_data.physical_get(idx_j, idx_i[0]);
				p[1] = i_data.physical_get(idx_j, idx_i[1]);
				p[2] = i_data.physical_get(idx_j, idx_i[2]);
				p[3] = i_data.physical_get(idx_j, idx_i[3]);
				p[4] = i_data.physical_get(idx_j, idx_i[4]);
				p[5] = i_data.physical_get(idx_j, idx_i[5]);

				double x0 = 0.;
				double x1 = 1.;
				double p0 = p[2];
				double p1 = p[3];

				double derx0  = (p[3] - p[1]) / 2.;
				double derx1  = (p[4] - p[2]) / 2.;
				double derxx0 = p[2] - 2. * p[1] + p[0];
				double derxx1 = p[5] - 2. * p[4] + p[3];

				double xx0 = x - x0;
				double xx0_2 = xx0 * xx0;
				double xx0_3 = xx0_2 * xx0;
				double xx1 = x - x1;
				double xx1_2 = xx1 * xx1;
				double xx1_3 = xx1_2 * xx1;
				double x1x0 = x1 -x0;
				double x1x0_2 = x1x0 * x1x0;
				double x1x0_3 = x1x0_2 * x1x0;
				double x1x0_4 = x1x0_3 * x1x0;
				double x1x0_5 = x1x0_4 * x1x0;
	
				q[kj] = p0 + derx0 * xx0 + .5 * derxx0 * xx0_2
				        + ( p1 - p0 - derx0 * x1x0 - .5 * derxx0 * x1x0_2) * xx0_3 / x1x0_3
				        + ( 3. * p0 - 3. * p1  + (2. * derx0 + derx1) * x1x0 + .5 * derxx0 * x1x0_2  ) * xx0_3 * xx1 / x1x0_4
				        + ( 6. * p1 - 6. * p0 - 3. * (derx0 + derx1) * x1x0 + .5 * (derxx1 - derxx0) * x1x0_2  ) * xx0_3 * xx1_2 / x1x0_5;

				idx_j = wrapPeriodic(idx_j+1, res[1]);
			}

			double y0 = 0.;
			double y1 = 1.;
			double q0 = q[2];
			double q1 = q[3];

			double dery0  = (q[3] - q[1]) / 2.;
			double dery1  = (q[4] - q[2]) / 2.;
			double deryy0 = q[2] - 2. * q[1] + q[0];
			double deryy1 = q[5] - 2. * q[4] + q[3];

			double yy0 = y - y0;
			double yy0_2 = yy0 * yy0;
			double yy0_3 = yy0_2 * yy0;
			double yy1 = y - y1;
			double yy1_2 = yy1 * yy1;
			double yy1_3 = yy1_2 * yy1;
			double y1y0 = y1 -y0;
			double y1y0_2 = y1y0 * y1y0;
			double y1y0_3 = y1y0_2 * y1y0;
			double y1y0_4 = y1y0_3 * y1y0;
			double y1y0_5 = y1y0_4 * y1y0;

			double value = q0 + dery0 * yy0 + .5 * deryy0 * yy0_2
			               + ( q1 - q0 - dery0 * y1y0 - .5 * deryy0 * y1y0_2) * yy0_3 / y1y0_3
			               + ( 3. * q0 - 3. * q1  + (2. * dery0 + dery1) * y1y0 + .5 * deryy0 * y1y0_2  ) * yy0_3 * yy1 / y1y0_4
			               + ( 6. * q1 - 6. * q0 - 3. * (dery0 + dery1) * y1y0 + .5 * (deryy1 - deryy0) * y1y0_2  ) * yy0_3 * yy1_2 / y1y0_5;

			o_data[pos_idx] = value;
		}
	}

public:
	void biquintic_scalar3(
			const PlaneData_Physical &i_data,			///< sampling data

			const ScalarDataArray &i_pos_x,		///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,		///< y positions of interpolation points

			double *o_data,						///< output values

			double i_shift_x = 0.0,				///< shift in x for staggered grids
			double i_shift_y = 0.0				///< shift in y for staggered grids
	)
	{
		assert(res[0] > 0);
		assert(cached_scale_factor[0] > 0);
		assert(i_pos_x.number_of_elements == i_pos_y.number_of_elements);

		std::size_t max_pos_idx = i_pos_x.number_of_elements;

#if SWEET_DEBUG
#if SWEET_THREADING_SPACE || SWEET_THREADING_TIME_REXI

		if (omp_get_num_threads() > 1)
			SWEETError("Are we already in parallel region? Threading race conditions likely!");
#endif
#endif

		// iterate over all positions in parallel
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t pos_idx = 0; pos_idx < max_pos_idx; pos_idx++)
		{
			/*
			 * load position to interpolate
			 * posx stores all x-coordinates of the arrival points
			 * posy stores all y-coordinates of the arrival points
			 *
			 * Both are arrays and matching array indices (pos_idx) below index the coordinates for the same point.
			 *
			 * Scale factor (Nx/dx, Ny/dy) maps from the physical space to the array space.
			 * The array space is from [0; N[
			 *
			 * shift_x/y is operating in array space. Hence, staggered grid can be
			 * implemented by setting this to 0.5 or -0.5
			 * for C grid, to interpolate given u data, use i_shift_x = 0.0,  i_shift_y = -0.5
			 *             to interpolate given v data, use i_shift_y = -0.5, i_shift_y = 0.0
			 *  pay attention to the negative shift, which is necessary because the staggered grids are positively shifted
			 *  and this shift has to be removed for the interpolation
			 */
			double pos_x = wrapPeriodic(i_pos_x.scalar_data[pos_idx]*cached_scale_factor[0] + i_shift_x, (double)res[0]);
			double pos_y = wrapPeriodic(i_pos_y.scalar_data[pos_idx]*cached_scale_factor[1] + i_shift_y, (double)res[1]);

			/**
			 * For the interpolation, we assume node-aligned values
			 *
			 * x0  x1  x2  x3
			 * |---|---|---|---
			 * 0   2   4   6    <- positions and associated values e.g. for domain size 8
			 */

			/**
			 * See http://www.paulinternet.nl/?page=bicubic
			 */
			// compute x/y position
			double x = pos_x - floor(pos_x);
			double y = pos_y - floor(pos_y);

			// precompute x-position indices since they are reused 4 times
			int idx_i[8];
			{

// TODO: Each 2nd wrapPeriodic is obsolete!!!
				int i = (int)pos_x-3;

				i = wrapPeriodic(i, res[0]);
				idx_i[0] = wrapPeriodic(i, res[0]);

				i = wrapPeriodic(i+1, res[0]);
				idx_i[1] = wrapPeriodic(i, res[0]);

				i = wrapPeriodic(i+1, res[0]);
				idx_i[2] = wrapPeriodic(i, res[0]);

				i = wrapPeriodic(i+1, res[0]);
				idx_i[3] = wrapPeriodic(i, res[0]);

				i = wrapPeriodic(i+1, res[0]);
				idx_i[4] = wrapPeriodic(i, res[0]);

				i = wrapPeriodic(i+1, res[0]);
				idx_i[5] = wrapPeriodic(i, res[0]);

				i = wrapPeriodic(i+1, res[0]);
				idx_i[6] = wrapPeriodic(i, res[0]);

				i = wrapPeriodic(i+1, res[0]);
				idx_i[7] = wrapPeriodic(i, res[0]);
			}

			/**
			 * iterate over rows and interpolate over the columns in the x direction
			 */
			// start at this row
			/////////////////int idx_j = wrapPeriodic((int)pos_y-2, res[1]);
			int idx_j = wrapPeriodic((int)pos_y - 3, res[1]);

			double q[8];
			for (int kj = 0; kj < 8; kj++)
			{
				double p[6];

				p[0] = i_data.physical_get(idx_j, idx_i[0]);
				p[1] = i_data.physical_get(idx_j, idx_i[1]);
				p[2] = i_data.physical_get(idx_j, idx_i[2]);
				p[3] = i_data.physical_get(idx_j, idx_i[3]);
				p[4] = i_data.physical_get(idx_j, idx_i[4]);
				p[5] = i_data.physical_get(idx_j, idx_i[5]);
				p[6] = i_data.physical_get(idx_j, idx_i[6]);
				p[7] = i_data.physical_get(idx_j, idx_i[7]);

				double x0 = 0.;
				double x1 = 1.;
				double p0 = p[3];
				double p1 = p[4];

				double derx0  = (p[4] - p[2]) / 2.;
				double derx1  = (p[5] - p[3]) / 2.;
				double derxx0 = 2. * p[3] - 5. * p[2] + 4. * p[1] - p[0];
				double derxx1 = 2. * p[4] - 5. * p[5] + 4. * p[6] - p[7];

				double xx0 = x - x0;
				double xx0_2 = xx0 * xx0;
				double xx0_3 = xx0_2 * xx0;
				double xx1 = x - x1;
				double xx1_2 = xx1 * xx1;
				double xx1_3 = xx1_2 * xx1;
				double x1x0 = x1 -x0;
				double x1x0_2 = x1x0 * x1x0;
				double x1x0_3 = x1x0_2 * x1x0;
				double x1x0_4 = x1x0_3 * x1x0;
				double x1x0_5 = x1x0_4 * x1x0;
	
				q[kj] = p0 + derx0 * xx0 + .5 * derxx0 * xx0_2
				        + ( p1 - p0 - derx0 * x1x0 - .5 * derxx0 * x1x0_2) * xx0_3 / x1x0_3
				        + ( 3. * p0 - 3. * p1  + (2. * derx0 + derx1) * x1x0 + .5 * derxx0 * x1x0_2  ) * xx0_3 * xx1 / x1x0_4
				        + ( 6. * p1 - 6. * p0 - 3. * (derx0 + derx1) * x1x0 + .5 * (derxx1 - derxx0) * x1x0_2  ) * xx0_3 * xx1_2 / x1x0_5;

				idx_j = wrapPeriodic(idx_j+1, res[1]);
			}

			double y0 = 0.;
			double y1 = 1.;
			double q0 = q[3];
			double q1 = q[4];

			double dery0  = (q[4] - q[2]) / 2.;
			double dery1  = (q[5] - q[3]) / 2.;
			double deryy0 = 2. * q[3] - 5. * q[2] + 4. * q[1] - q[0];
			double deryy1 = 2. * q[4] - 5. * q[5] + 4. * q[6] - q[7];

			double yy0 = y - y0;
			double yy0_2 = yy0 * yy0;
			double yy0_3 = yy0_2 * yy0;
			double yy1 = y - y1;
			double yy1_2 = yy1 * yy1;
			double yy1_3 = yy1_2 * yy1;
			double y1y0 = y1 -y0;
			double y1y0_2 = y1y0 * y1y0;
			double y1y0_3 = y1y0_2 * y1y0;
			double y1y0_4 = y1y0_3 * y1y0;
			double y1y0_5 = y1y0_4 * y1y0;

			double value = q0 + dery0 * yy0 + .5 * deryy0 * yy0_2
			               + ( q1 - q0 - dery0 * y1y0 - .5 * deryy0 * y1y0_2) * yy0_3 / y1y0_3
			               + ( 3. * q0 - 3. * q1  + (2. * dery0 + dery1) * y1y0 + .5 * deryy0 * y1y0_2  ) * yy0_3 * yy1 / y1y0_4
			               + ( 6. * q1 - 6. * q0 - 3. * (dery0 + dery1) * y1y0 + .5 * (deryy1 - deryy0) * y1y0_2  ) * yy0_3 * yy1_2 / y1y0_5;

			o_data[pos_idx] = value;
		}
	}



public:
	const PlaneData_Physical biquintic_scalar(
			PlaneData_Physical &i_data,				///< sampling data

			const ScalarDataArray &i_pos_x,				///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,				///< y positions of interpolation points
			//PlaneData* i_pos[2],	///< sampling position
			double i_shift_x = 0.0,
			double i_shift_y = 0.0
	)
	{
		PlaneData_Physical out(planeDataConfig);
		biquintic_scalar(i_data, i_pos_x, i_pos_y, out, i_shift_x, i_shift_y);
		return out;
	}

public:
	void biquintic_scalar(
			const PlaneData_Physical &i_data,			///< sampling data

			const ScalarDataArray &i_pos_x,		///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,		///< y positions of interpolation points

			PlaneData_Physical &o_data,					///< output values

			double i_shift_x = 0.0,				///< shift in x for staggered grids
			double i_shift_y = 0.0				///< shift in y for staggered grids
	)
	{
		assert(res[0] > 0);
		assert(cached_scale_factor[0] > 0);
		assert(i_pos_x.number_of_elements == i_pos_y.number_of_elements);
		assert(i_pos_x.number_of_elements == o_data.planeDataConfig->physical_array_data_number_of_elements);

		biquintic_scalar(i_data, i_pos_x, i_pos_y, o_data.physical_space_data, i_shift_x, i_shift_y);


	}


public:
	const PlaneData_Physical bi_interp_scalar(
			PlaneData_Physical &i_data,				///< sampling data

			const ScalarDataArray &i_pos_x,				///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,				///< y positions of interpolation points

			int order,

			//PlaneData* i_pos[2],	///< sampling position
			double i_shift_x = 0.0,
			double i_shift_y = 0.0

	)
	{
		PlaneData_Physical out(planeDataConfig);
		bi_interp_scalar(i_data, i_pos_x, i_pos_y, out, order, i_shift_x, i_shift_y);
		return out;
	}

public:
	PlaneData_Spectral bi_interp_scalar(
			const PlaneData_Spectral &i_data,			///< sampling data

			const ScalarDataArray &i_pos_x,		///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,		///< y positions of interpolation points

			int order,

			double i_shift_x = 0.0,				///< shift in x for staggered grids
			double i_shift_y = 0.0				///< shift in y for staggered grids

	)
	{

		PlaneData_Physical i_data_phys = i_data.toPhys();
		PlaneData_Spectral out(i_data.planeDataConfig);

		PlaneData_Physical out_phys = bi_interp_scalar(i_data_phys, i_pos_x, i_pos_y, order, i_shift_x, i_shift_y);

		out.loadPlaneDataPhysical(out_phys);

		return out;

	}



public:
	void bi_interp_scalar(
			const PlaneData_Physical &i_data,			///< sampling data

			const ScalarDataArray &i_pos_x,		///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,		///< y positions of interpolation points

			PlaneData_Physical &o_data,					///< output values

			int order,

			double i_shift_x = 0.0,				///< shift in x for staggered grids
			double i_shift_y = 0.0				///< shift in y for staggered grids

	)
	{
		assert(res[0] > 0);
		assert(cached_scale_factor[0] > 0);
		assert(i_pos_x.number_of_elements == i_pos_y.number_of_elements);
		assert(i_pos_x.number_of_elements == o_data.planeDataConfig->physical_array_data_number_of_elements);

                if (order == 1)
			bilinear_scalar(i_data, i_pos_x, i_pos_y, o_data.physical_space_data, i_shift_x, i_shift_y);
		else if (order == 32)
			bicubic_scalar_NEW2(i_data, i_pos_x, i_pos_y, o_data.physical_space_data, i_shift_x, i_shift_y);
		else if (order == 34)
			bicubic_scalar_NEW4(i_data, i_pos_x, i_pos_y, o_data.physical_space_data, i_shift_x, i_shift_y);
		else if (order == 36)
			bicubic_scalar_NEW6(i_data, i_pos_x, i_pos_y, o_data.physical_space_data, i_shift_x, i_shift_y);
		else if (order == 38)
			bicubic_scalar_NEW8(i_data, i_pos_x, i_pos_y, o_data.physical_space_data, i_shift_x, i_shift_y);
		else if (order == -3)
			bicubic_scalar_Lagrangian(i_data, i_pos_x, i_pos_y, o_data.physical_space_data, i_shift_x, i_shift_y);
		else if (order == 5)
			biquintic_scalar(i_data, i_pos_x, i_pos_y, o_data.physical_space_data, i_shift_x, i_shift_y);
		else if (order == 55)
			biquintic_scalar2(i_data, i_pos_x, i_pos_y, o_data.physical_space_data, i_shift_x, i_shift_y);
		else if (order == 555)
			biquintic_scalar3(i_data, i_pos_x, i_pos_y, o_data.physical_space_data, i_shift_x, i_shift_y);
		else if (order == -5)
			biquintic_scalar_Lagrangian(i_data, i_pos_x, i_pos_y, o_data.physical_space_data, i_shift_x, i_shift_y);
		else if (order == -7)
			biseptic_scalar_Lagrangian(i_data, i_pos_x, i_pos_y, o_data.physical_space_data, i_shift_x, i_shift_y);
		else
			SWEETError("This interpolation order is not implemented");

	}


};




#endif /* SRC_INCLUDE_SWEET_PLANEDATASAMPLER_HPP_ */
