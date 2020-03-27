/*
 * SWEETMath.hpp
 *
 *  Created on: Mar 27, 2020
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */

#include <sweet/ScalarDataArray.hpp>

#ifndef SRC_INCLUDE_SWEET_MATH_HPP_
#define SRC_INCLUDE_SWEET_MATH_HPP_



class SWEETMath
{
public:
	static
	void normalize_vec3(
		const double i_vec[3],
		double o_vec[3]
	)
	{
		double s = 1.0/std::sqrt(i_vec[0]*i_vec[0] + i_vec[1]*i_vec[1] + i_vec[2]*i_vec[2]);

		o_vec[0] = i_vec[0] * s;
		o_vec[1] = i_vec[1] * s;
		o_vec[2] = i_vec[2] * s;
	}



	static
	void normalize_vec3(
		double i_vec[3]
	)
	{
		double s = 1.0/std::sqrt(i_vec[0]*i_vec[0] + i_vec[1]*i_vec[1] + i_vec[2]*i_vec[2]);

		i_vec[0] *= s;
		i_vec[1] *= s;
		i_vec[2] *= s;
	}


	static
	void normalize(
			ScalarDataArray &io_v0,
			ScalarDataArray &io_v1,
			ScalarDataArray &io_v2
	)
	{
		ScalarDataArray norm = (io_v0*io_v0 + io_v1*io_v1 + io_v2*io_v2).inv_sqrt();

		io_v0 *= norm;
		io_v1 *= norm;
		io_v2 *= norm;

	}


	/*
	 * Normalize by using threshold
	 *
	 * If the L2^2 value is smaller than the threshold, the norm is set to 1 instead of infinity
	 */
	static
	void normalize_threshold(
			ScalarDataArray &io_v0,
			ScalarDataArray &io_v1,
			ScalarDataArray &io_v2,
			double i_threshold
	)
	{
		double threshold2 = i_threshold*i_threshold;

		ScalarDataArray length2 = io_v0*io_v0 + io_v1*io_v1 + io_v2*io_v2;

		length2.update_lambda_array_indices(
				[&](int i, double &io_data)
				{
					if (io_data < threshold2)
						io_data = 1.0;
				}
		);
		ScalarDataArray norm = length2.inv_sqrt();

		io_v0 *= norm;
		io_v1 *= norm;
		io_v2 *= norm;

	}



	static
	ScalarDataArray length(
			const ScalarDataArray &io_v0,
			const ScalarDataArray &io_v1,
			const ScalarDataArray &io_v2
	)
	{
		return (io_v0*io_v0 + io_v1*io_v1 + io_v2*io_v2).sqrt();
	}



	static
	void cross_prod_v3(double w[3], double v[3], double ret[3])
	{
		ret[0] = w[1]*v[2] - w[2]*v[1];
		ret[1] = w[2]*v[0] - w[0]*v[2];
		ret[2] = w[0]*v[1] - w[1]*v[0];
	}


	static
	void cross_prod(
			const ScalarDataArray &w0,
			const ScalarDataArray &w1,
			const ScalarDataArray &w2,

			const ScalarDataArray &v0,
			const ScalarDataArray &v1,
			const ScalarDataArray &v2,

			ScalarDataArray &o_ret0,
			ScalarDataArray &o_ret1,
			ScalarDataArray &o_ret2
	)
	{
		o_ret0 = w1*v2 - w2*v1;
		o_ret1 = w2*v0 - w0*v2;
		o_ret2 = w0*v1 - w1*v0;
	}


	inline
	static
	void latlon_to_cartesian(
			const ScalarDataArray &i_lon,
			const ScalarDataArray &i_lat,
			ScalarDataArray &o_x,
			ScalarDataArray &o_y,
			ScalarDataArray &o_z
	)
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t i = 0; i < i_lon.number_of_elements; i++)
		{
			o_x.scalar_data[i] = std::cos(i_lon.scalar_data[i])*std::cos(i_lat.scalar_data[i]);
			o_y.scalar_data[i] = std::sin(i_lon.scalar_data[i])*std::cos(i_lat.scalar_data[i]);
			o_z.scalar_data[i] = std::sin(i_lat.scalar_data[i]);
		}
	}



	inline
	static
	void latlon_to_cartesian(
			const double i_lon,
			const double &i_lat,
			double o_ret[3]
	)
	{
		o_ret[0] = std::cos(i_lon)*std::cos(i_lat);
		o_ret[1] = std::sin(i_lon)*std::cos(i_lat);
		o_ret[2] = std::sin(i_lat);
	}



	inline
	static
	void latlon_velocity_to_cartesian_velocity(
			const ScalarDataArray &i_lon,
			const ScalarDataArray &i_lat,
			const ScalarDataArray &i_vel_lon,
			const ScalarDataArray &i_vel_lat,
			ScalarDataArray *o_v_x,
			ScalarDataArray *o_v_y,
			ScalarDataArray *o_v_z
	)
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t i = 0; i < i_lon.number_of_elements; i++)
		{
			o_v_x->scalar_data[i] = -i_vel_lon.scalar_data[i]*std::sin(i_lon.scalar_data[i]) - i_vel_lat.scalar_data[i]*std::cos(i_lon.scalar_data[i])*std::sin(i_lat.scalar_data[i]);
			o_v_y->scalar_data[i] = i_vel_lon.scalar_data[i]*std::cos(i_lon.scalar_data[i]) - i_vel_lat.scalar_data[i]*std::sin(i_lon.scalar_data[i])*std::sin(i_lat.scalar_data[i]);
			o_v_z->scalar_data[i] = i_vel_lat.scalar_data[i]*std::cos(i_lat.scalar_data[i]);
		}
	}



	static
	void cartesian_to_latlon(
			const ScalarDataArray &i_x,
			const ScalarDataArray &i_y,
			const ScalarDataArray &i_z,
			ScalarDataArray &o_lon,
			ScalarDataArray &o_lat
	)
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t i = 0; i < i_x.number_of_elements; i++)
		{
			/*
			 * Make sure that coordinates are in valid range
			 */
#if 0
			i_x.scalar_data[i] = std::min(1., i_x.scalar_data[i]);
			i_x.scalar_data[i] = std::max(-1., i_x.scalar_data[i]);

			i_y.scalar_data[i] = std::min(1., i_y.scalar_data[i]);
			i_y.scalar_data[i] = std::max(-1., i_y.scalar_data[i]);

			i_z.scalar_data[i] = std::min(1., i_z.scalar_data[i]);
			i_z.scalar_data[i] = std::max(-1., i_z.scalar_data[i]);
#endif

			/*
			 * Now compute the angles
			 */
			o_lon.scalar_data[i] = std::atan(i_y.scalar_data[i]/i_x.scalar_data[i]);

#if SWEET_DEBUG
			if (	(std::isnan(o_lon.scalar_data[i]) != 0) ||
					(std::abs(std::isinf(o_lon.scalar_data[i])) == 1)
			)
			{
				std::cout << "Found nan/inf at position " << i << std::endl;
				std::cout << "o_lon value: " << o_lon.scalar_data[i] << std::endl;
				std::cout << "   atan(" << i_y.scalar_data[i] << ", " << i_x.scalar_data[i] << ")" << std::endl;
				FatalError("EXIT");
			}
#endif

			if (i_x.scalar_data[i] < 0)
				o_lon.scalar_data[i] += M_PI;
			else if (i_y.scalar_data[i] < 0)
				o_lon.scalar_data[i] += M_PI*2.0;

			o_lat.scalar_data[i] = std::acos(-i_z.scalar_data[i]) - M_PI*0.5;

#if SWEET_DEBUG
			if (
					(std::isnan(o_lat.scalar_data[i]) != 0) ||
					(std::abs(std::isinf(o_lat.scalar_data[i])) == 1)
			)
			{
				std::cout << "Found nan/inf at position " << i << std::endl;
				std::cout << "o_lat value: " << o_lat.scalar_data[i] << std::endl;
				std::cout << "   acos(" << i_z.scalar_data[i] << ")" << std::endl;
				FatalError("EXIT");
			}
#endif
		}
	}



	/**
	 * Rotate a point around an axis for the given rotation angle.
	 */
	static
	void rotate_3d(
		double i_pos_start[3],
		double i_rotation_angle,
		double i_rotation_axis[3],
		double i_pos_new[3]
	)
	{
		/*
		 * Based on source code from website:
		 * http://paulbourke.net/geometry/rotate/
		 * http://paulbourke.net/geometry/rotate/source.c (by Ronald Goldman)
		 */

	   normalize_vec3(i_rotation_axis);
	   double costheta = std::cos(i_rotation_angle);
	   double sintheta = std::sin(i_rotation_angle);

	   i_pos_new[0] = (costheta + (1 - costheta) * i_rotation_axis[0] * i_rotation_axis[0]) * i_pos_start[0];
	   i_pos_new[0] += ((1 - costheta) * i_rotation_axis[0] * i_rotation_axis[1] - i_rotation_axis[2] * sintheta) * i_pos_start[1];
	   i_pos_new[0] += ((1 - costheta) * i_rotation_axis[0] * i_rotation_axis[2] + i_rotation_axis[1] * sintheta) * i_pos_start[2];

	   i_pos_new[1] = ((1 - costheta) * i_rotation_axis[0] * i_rotation_axis[1] + i_rotation_axis[2] * sintheta) * i_pos_start[0];
	   i_pos_new[1] += (costheta + (1 - costheta) * i_rotation_axis[1] * i_rotation_axis[1]) * i_pos_start[1];
	   i_pos_new[1] += ((1 - costheta) * i_rotation_axis[1] * i_rotation_axis[2] - i_rotation_axis[0] * sintheta) * i_pos_start[2];

	   i_pos_new[2] = ((1 - costheta) * i_rotation_axis[0] * i_rotation_axis[2] - i_rotation_axis[1] * sintheta) * i_pos_start[0];
	   i_pos_new[2] += ((1 - costheta) * i_rotation_axis[1] * i_rotation_axis[2] + i_rotation_axis[0] * sintheta) * i_pos_start[1];
	   i_pos_new[2] += (costheta + (1 - costheta) * i_rotation_axis[2] * i_rotation_axis[2]) * i_pos_start[2];
	}


	static
	ScalarDataArray cos(const ScalarDataArray &i_values)
	{
		ScalarDataArray retval = i_values;
		retval.update_lambda_array_indices(
			[](int, double &io_data)
			{
				io_data = std::cos(io_data);
			}
		);

		return retval;
	}



	static
	ScalarDataArray sin(const ScalarDataArray &i_values)
	{
		ScalarDataArray retval = i_values;
		retval.update_lambda_array_indices(
			[](int, double &io_data)
			{
				io_data = std::sin(io_data);
			}
		);

		return retval;
	}


	/**
	 * Rotate a point around an axis for the given rotation angle.
	 */
	static
	void rotate_3d(
			const ScalarDataArray &i_pos_start_0,
			const ScalarDataArray &i_pos_start_1,
			const ScalarDataArray &i_pos_start_2,
			const ScalarDataArray &i_rotation_angle,
			const ScalarDataArray &i_rotation_axis_0,
			const ScalarDataArray &i_rotation_axis_1,
			const ScalarDataArray &i_rotation_axis_2,
			ScalarDataArray &o_pos_new_0,
			ScalarDataArray &o_pos_new_1,
			ScalarDataArray &o_pos_new_2
	)
	{
		/*
		 * Based on source code from website:
		 * http://paulbourke.net/geometry/rotate/
		 * http://paulbourke.net/geometry/rotate/source.c (by Ronald Goldman)
		 */

		ScalarDataArray rotation_axis_0 = i_rotation_axis_0;
		ScalarDataArray rotation_axis_1 = i_rotation_axis_1;
		ScalarDataArray rotation_axis_2 = i_rotation_axis_2;

		normalize(
				rotation_axis_0,
				rotation_axis_1,
				rotation_axis_2
			);

		ScalarDataArray cos_theta = cos(i_rotation_angle);
		ScalarDataArray sin_theta = sin(i_rotation_angle);

		o_pos_new_0 = (cos_theta + (1 - cos_theta) * rotation_axis_0 * rotation_axis_0) * i_pos_start_0;
		o_pos_new_0 += ((1 - cos_theta) * rotation_axis_0 * rotation_axis_1 - rotation_axis_2 * sin_theta) * i_pos_start_1;
		o_pos_new_0 += ((1 - cos_theta) * rotation_axis_0 * rotation_axis_2 + rotation_axis_1 * sin_theta) * i_pos_start_2;

		o_pos_new_1 = ((1 - cos_theta) * rotation_axis_0 * rotation_axis_1 + rotation_axis_2 * sin_theta) * i_pos_start_0;
		o_pos_new_1 += (cos_theta + (1 - cos_theta) * rotation_axis_1 * rotation_axis_1) * i_pos_start_1;
		o_pos_new_1 += ((1 - cos_theta) * rotation_axis_1 * rotation_axis_2 - rotation_axis_0 * sin_theta) * i_pos_start_2;

		o_pos_new_2 = ((1 - cos_theta) * rotation_axis_0 * rotation_axis_2 - rotation_axis_1 * sin_theta) * i_pos_start_0;
		o_pos_new_2 += ((1 - cos_theta) * rotation_axis_1 * rotation_axis_2 + rotation_axis_0 * sin_theta) * i_pos_start_1;
		o_pos_new_2 += (cos_theta + (1 - cos_theta) * rotation_axis_2 * rotation_axis_2) * i_pos_start_2;
	}




	/**
	 * Rotate a point around an axis for the given rotation angle.
	 *
	 * The rotation axis is assumed to be already normalized
	 */
	static
	void rotate_3d_normalized_rotation_axis(
		double i_pos_start[3],
		double i_rotation_angle,
		double i_rotation_axis_normalized[3],
		double i_pos_new[3]
	)
	{
		/*
		 * Based on source code from website:
		 * http://paulbourke.net/geometry/rotate/
		 * http://paulbourke.net/geometry/rotate/source.c (by Ronald Goldman)
		 */

	   double costheta = std::cos(i_rotation_angle);
	   double sintheta = std::sin(i_rotation_angle);

	   i_pos_new[0] = (costheta + (1 - costheta) * i_rotation_axis_normalized[0] * i_rotation_axis_normalized[0]) * i_pos_start[0];
	   i_pos_new[0] += ((1 - costheta) * i_rotation_axis_normalized[0] * i_rotation_axis_normalized[1] - i_rotation_axis_normalized[2] * sintheta) * i_pos_start[1];
	   i_pos_new[0] += ((1 - costheta) * i_rotation_axis_normalized[0] * i_rotation_axis_normalized[2] + i_rotation_axis_normalized[1] * sintheta) * i_pos_start[2];

	   i_pos_new[1] = ((1 - costheta) * i_rotation_axis_normalized[0] * i_rotation_axis_normalized[1] + i_rotation_axis_normalized[2] * sintheta) * i_pos_start[0];
	   i_pos_new[1] += (costheta + (1 - costheta) * i_rotation_axis_normalized[1] * i_rotation_axis_normalized[1]) * i_pos_start[1];
	   i_pos_new[1] += ((1 - costheta) * i_rotation_axis_normalized[1] * i_rotation_axis_normalized[2] - i_rotation_axis_normalized[0] * sintheta) * i_pos_start[2];

	   i_pos_new[2] = ((1 - costheta) * i_rotation_axis_normalized[0] * i_rotation_axis_normalized[2] - i_rotation_axis_normalized[1] * sintheta) * i_pos_start[0];
	   i_pos_new[2] += ((1 - costheta) * i_rotation_axis_normalized[1] * i_rotation_axis_normalized[2] + i_rotation_axis_normalized[0] * sintheta) * i_pos_start[1];
	   i_pos_new[2] += (costheta + (1 - costheta) * i_rotation_axis_normalized[2] * i_rotation_axis_normalized[2]) * i_pos_start[2];
	}


	/**
	 * Rotate a point around an axis for the given rotation angle.
	 */
	static
	void rotate_3d_normalized_rotation_axis(
			const ScalarDataArray &i_pos_start_0,
			const ScalarDataArray &i_pos_start_1,
			const ScalarDataArray &i_pos_start_2,
			const ScalarDataArray &i_rotation_angle,
			const ScalarDataArray &i_rotation_axis_0,
			const ScalarDataArray &i_rotation_axis_1,
			const ScalarDataArray &i_rotation_axis_2,
			ScalarDataArray &o_pos_new_0,
			ScalarDataArray &o_pos_new_1,
			ScalarDataArray &o_pos_new_2
	)
	{
		/*
		 * Based on source code from website:
		 * http://paulbourke.net/geometry/rotate/
		 * http://paulbourke.net/geometry/rotate/source.c (by Ronald Goldman)
		 */

		ScalarDataArray cos_theta = cos(i_rotation_angle);
		ScalarDataArray sin_theta = sin(i_rotation_angle);

		o_pos_new_0 = (cos_theta + (1 - cos_theta) * i_rotation_axis_0 * i_rotation_axis_0) * i_pos_start_0;
		o_pos_new_0 += ((1 - cos_theta) * i_rotation_axis_0 * i_rotation_axis_1 - i_rotation_axis_2 * sin_theta) * i_pos_start_1;
		o_pos_new_0 += ((1 - cos_theta) * i_rotation_axis_0 * i_rotation_axis_2 + i_rotation_axis_1 * sin_theta) * i_pos_start_2;

		o_pos_new_1 = ((1 - cos_theta) * i_rotation_axis_0 * i_rotation_axis_1 + i_rotation_axis_2 * sin_theta) * i_pos_start_0;
		o_pos_new_1 += (cos_theta + (1 - cos_theta) * i_rotation_axis_1 * i_rotation_axis_1) * i_pos_start_1;
		o_pos_new_1 += ((1 - cos_theta) * i_rotation_axis_1 * i_rotation_axis_2 - i_rotation_axis_0 * sin_theta) * i_pos_start_2;

		o_pos_new_2 = ((1 - cos_theta) * i_rotation_axis_0 * i_rotation_axis_2 - i_rotation_axis_1 * sin_theta) * i_pos_start_0;
		o_pos_new_2 += ((1 - cos_theta) * i_rotation_axis_1 * i_rotation_axis_2 + i_rotation_axis_0 * sin_theta) * i_pos_start_1;
		o_pos_new_2 += (cos_theta + (1 - cos_theta) * i_rotation_axis_2 * i_rotation_axis_2) * i_pos_start_2;
	}


};


#endif
