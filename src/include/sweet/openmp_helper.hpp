/*
 * openmp_helper.hpp
 *
 *  Created on: 20 Jul 2015
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */
#ifndef SRC_EXAMPLES_OPENMP_HELPER_HPP_
#define SRC_EXAMPLES_OPENMP_HELPER_HPP_



/**
 * This is a class to overcome SIMD instruction issues with older GNU compilers
 */

#define OMP_SCHEDULE	schedule(static)

#if !SWEET_THREADING

	#define SWEET_OMP_PARALLEL_FOR_SIMD
	#define SWEET_OMP_PARALLEL_FOR

	#define SWEET_THREADING_SPACE_PARALLEL_FOR
	#define SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
	#define SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2

	#define SWEET_THREADING_TIME_PARALLEL_FOR
	#define SWEET_THREADING_TIME_PARALLEL_FOR_SIMD
	#define SWEET_THREADING_TIME_PARALLEL_FOR_SIMD_COLLAPSE2

	#define PROC_BIND_CLOSE

#else

	#define SWEET_OMP_PARALLEL_FOR _Pragma("omp parallel for")

	// http://gcc.gnu.org/onlinedocs/cpp/Common-Predefined-Macros.html
	#if __GNUC__ >= 7
		#define PROC_BIND_CLOSE proc_bind(close)
	#else
		#define PROC_BIND_CLOSE 
	#endif

	#if SWEET_SIMD_ENABLE
		#define SWEET_OMP_PARALLEL_FOR_SIMD _Pragma("omp parallel for simd")
		#define SWEET_OMP_PARALLEL_FOR_SIMD_COLLAPSE2 _Pragma("omp parallel for simd collapse(2)")
	#else
		#define SWEET_OMP_PARALLEL_FOR_SIMD _Pragma("omp parallel for")
	#endif

	#if SWEET_THREADING_SPACE
		#define SWEET_THREADING_SPACE_PARALLEL_FOR SWEET_OMP_PARALLEL_FOR
		#define SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD SWEET_OMP_PARALLEL_FOR_SIMD
		#define SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2 SWEET_OMP_PARALLEL_FOR_SIMD_COLLAPSE2
	#else
		#define SWEET_THREADING_SPACE_PARALLEL_FOR
		#define SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		#define SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2
	#endif


	#if SWEET_THREADING_TIME
		#define SWEET_THREADING_TIME_PARALLEL_FOR_SIMD SWEET_OMP_PARALLEL_FOR_SIMD
		#define SWEET_THREADING_TIME_PARALLEL_FOR SWEET_OMP_PARALLEL_FOR
		#define SWEET_THREADING_TIME_PARALLEL_FOR_SIMD_COLLAPSE2 SWEET_OMP_PARALLEL_FOR_SIMD_COLLAPSE2
	#else
		#define SWEET_THREADING_TIME_PARALLEL_FOR_SIMD
		#define SWEET_THREADING_TIME_PARALLEL_FOR
		#define SWEET_THREADING_TIME_PARALLEL_FOR_SIMD_COLLAPSE2 SWEET_OMP_PARALLEL_FOR_SIMD_COLLAPSE2
	#endif


#endif



#endif /* SRC_EXAMPLES_OPENMP_HELPER_HPP_ */
