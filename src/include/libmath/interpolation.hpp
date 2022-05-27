/*
 * interpolation.cpp
 *
 *  Created on: 18 March 2020
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */


#ifndef SWEET_INTERPOLATE_HPP
#define SWEET_INTERPOLATE_HPP

/**
 * Compute Lagrange interpolation for given interpolation points and interpolant values
 * 
 * https://de.wikipedia.org/wiki/Polynominterpolation#Lagrangesche_Interpolationsformel
 *
 * Here we need to provide the nonequidistantly coordinates of the interpolation points.
 */
template <int N>
double interpolation_lagrange_nonequidistant(
    double *x,    /// interpolation points
    double *y,    /// interpolation values
    double x_sample /// sample position
)
{
    double retval = 0;
    for (int i = 0; i < N; i++)
    {
        double denom = 1;
        double nom = 1;

        for (int j = 0; j < N; j++)
        {
            if (i == j)
                continue;

            nom *= x_sample - x[j];
            denom *= x[i] - x[j];
        }

        retval += y[i]*nom/denom;
    }

    return retval;
}




/**
 * Compute Lagrange interpolation for given interpolation points and interpolant values
 *
 * https://de.wikipedia.org/wiki/Polynominterpolation#Lagrangesche_Interpolationsformel
 *
 * Here we assume the equidistantly spaced points starting at 0
 */
template <int N>
double interpolation_lagrange_equidistant(
    double *y,    /// interpolation values
    double x_sample /// sample position
)
{
    double retval = 0;
    for (int i = 0; i < N; i++)
    {
        double denom = 1;
        double nom = 1;

        for (int j = 0; j < N; j++)
        {
            if (i == j)
                continue;

            nom *= x_sample - (double)j;
            denom *= (double)(i - j);
        }

        retval += y[i]*nom/denom;
    }

    return retval;
}



inline double factorial(int n)
{
    double p = 1.;
    for (int i = 2; i <= n; i++)
        p *= i;
    return p;
}

/**
 * Derivative approximation to be used in Hermite interpolation
 */
inline double derivative(
	int n,      // derivative
	int p,      // derivative order
	std::vector<double> x,  // points
	std::vector<double> y  // values
)
{

	double d;

	if (n == 1)
	{
		if (p == 2)
		{
			d = (y[2] - y[0]) / (x[2] - x[0]);
		}
	}
	else
		SWEETError("Not implemented! n = " + std::to_string(n) + ", p = " + std::to_string(p));

	return d;
/////	std::vector<double> c = {}; // coefficients
/////
/////	if (n == 1)
/////	{
/////		if (p == 2)
/////			c = 	{
/////					-.5,
/////					0.,
/////					.5
/////				};
/////		else if (p == 4)
/////			c = 	{
/////					1./12.,
/////					-2./3.,
/////					0.,
/////					2./3.,
/////					-1./12.
/////				}
/////		else if (p == 6)
/////			c =	{
/////					-1./60.,
/////					3./20.,
/////					-3./4.,
/////					0.,
/////					3./4.,
/////					-3./20.,
/////					1./60.
/////				}
/////		else if (p == 8)
/////			c =	{
/////					1./280.,
/////					-4./105.,
/////					1./5.,
/////					-4./5.,
/////					0.,
/////					4./5.,
/////					-1./5,
/////					4./105.,
/////					-1./280.
/////				}
/////	}
/////	else if (n == 2)
/////	{
/////		if (p == 2)
/////		{
/////			c = 	{
/////					1.,
/////					-2.,
/////					1.
/////				}
/////		}
/////		else if (p == 4)
/////		{
/////			c =	{
/////					-1./12.,
/////					4./3.,
/////					-5./2.,
/////					4./3.,
/////					-1./12.
/////				}
/////		}
/////	}




}


/**
 * Differential operator for computing recursively
 * May be too expensive; better to compute coefficients only once
 *
 */
inline double diff_hermite(
	std::vector<double> x,		// entries for computing differences
	std::vector<double> y,		// entries for computing differences
        std::vector<double> x_all,	// (all) entries for computing derivatives
        std::vector<double> y_all,	// (all) entries for computing derivatives
	int n,
	int p,
	int N
)
{

	//////std::cout << std::endl;
	//////std::cout << "X" << std::endl;
	//////for (size_t i = 0; i < x.size(); i++)
	//////	std::cout << x[i] << " ";
	//////std::cout << std::endl;
	//////std::cout << "Y" << std::endl;
	//////for (size_t i = 0; i < y.size(); i++)
	//////	std::cout << y[i] << " ";
	//////std::cout << std::endl;

	// if only one element, return it
	if (y.size() == 1)
		return y[0];

	bool all_equal = true;
	for (size_t i = 1; i < x.size(); i++)
		if (x[i] != x[0])
		{
			all_equal = false;
			break;
		}

	// if all xs are the same, return the derivative
	//if ( std::equal(x.begin() + 1, x.end(), x.begin()) )
	if ( all_equal )
	{
		///std::cout << "ALL EQUAL" << std::endl;
		//assert (n == (int)y.size() + 1);

		// get derivative order
		int der_order = (int)y.size() - 1;

		// get stencil for computing the derivative
		// first: check if remaining point is x0 or x1
		int b = (int)(N/2) - 1;
		std::vector<double> x_stencil;
		std::vector<double> y_stencil;
		if (x[0] == x_all[b]) // x0
		{
			x_stencil = std::vector<double>(x_all.begin(), x_all.end() - 1); // check this!
			y_stencil = std::vector<double>(y_all.begin(), y_all.end() - 1); // check this!
		}
		else if (x[0] == x_all[b+1]) // x1
		{
			x_stencil = std::vector<double>(x_all.begin() + 1, x_all.end()); // check this!
			y_stencil = std::vector<double>(y_all.begin() + 1, y_all.end()); // check this!
		}
		else
			SWEETError("Wrong x " + std::to_string(x[0]) + " " + std::to_string(x_all[b]) + " " + std::to_string(x_all[b + 1]));

		// compute the derivative
		///std::cout << "COMPUTING DERIVATIVE with " << std::endl;
		///std::cout << "SIZE VECTORS : " << x_stencil.size() << " " << y_stencil.size() << std::endl;
		///std::cout << "x " << x_stencil[0] << " " << x_stencil[1] << " " << x_stencil[2] << std::endl;
		///std::cout << "y " << y_stencil[0] << " " << y_stencil[1] << " " << y_stencil[2] << std::endl;
		///std::cout << "DERIVATIVE:  " << der_order << " " << derivative(der_order, p, x_stencil, y_stencil) << " " << factorial(der_order) << std::endl;
		return derivative(der_order, p, x_stencil, y_stencil) / factorial(der_order);
	}

	// else, recursive difference
	///std::cout << "RECURSIVE" << std::endl;
	std::vector<double> x0 = std::vector<double>(x.begin(), x.end() - 1);
	std::vector<double> y0 = std::vector<double>(y.begin(), y.end() - 1);
	std::vector<double> x1 = std::vector<double>(x.begin() + 1, x.end());
	std::vector<double> y1 = std::vector<double>(y.begin() + 1, y.end());
	///std::cout << "DIFFs " << diff_hermite(x1, y1, x_all, y_all, n, p, N) << " " << diff_hermite(x0, y0, x_all, y_all, n, p, N) << " " << x.back() - x[0]  << std::endl;
	return ( diff_hermite(x1, y1, x_all, y_all, n, p, N) - diff_hermite(x0, y0, x_all, y_all, n, p, N) ) / ( x.back() - x[0] );

}



/**
 * Compute n-th Hermite interpolation for two given interpolation points
 * on which we know the interpolation value and the first (n-1)/2 derivatives
 * A given p-th order centered approximation to the derivatives is considered.
 * A total of N = p + 2 + 2 * ( int(n-2)/2 )
 *
 * https://de.wikipedia.org/wiki/Polynominterpolation#Lagrangesche_Interpolationsformel
 *
 * Here we need to provide the nonequidistantly coordinates of the interpolation points.
 */
template <int N, int n, int p>
double interpolation_hermite(
    double *y,       /// interpolation values
    double x_sample, /// sample position
    bool equidistant,
    double *x = nullptr
)
{

    ///std::cout << std::endl << std::endl << "Starting" << std::endl;
    ///std::cout << "EQUIDISTANT " << equidistant << std::endl;

    // only odd order polynomials!
    assert (n%2 == 1);
    assert (N == p + 2 + 2 * ( (n - 2) / 2 ));
    assert (sizeof(y) / sizeof(y[0]) == N);

    // convert x and y to vectors
    std::vector<double> x_vec = {};
    std::vector<double> y_vec = {};
    for (size_t i = 0; i < N; i++)
    {
        if (equidistant)
            x_vec.push_back(i);
        else
            x_vec.push_back(x[i]);
        y_vec.push_back(y[i]);
    }

    /////std::cout << "x_vec" << std::endl;
    /////for (size_t i = 0; i < x_vec.size(); i++)
    /////    std::cout << x_vec[i] << " ";
    /////std::cout << std::endl;
    /////std::cout << "y_vec" << std::endl;
    /////for (size_t i = 0; i < y_vec.size(); i++)
    /////    std::cout << y_vec[i] << " ";
    /////std::cout << std::endl;


    // vector with repeated entries indicating derivatives
    // We repeat the points x0 and x1
    // located at position N/2-1 and N/2
    std::vector<double> wx = {};
    std::vector<double> wy = {};
    int b = (int)(N/2) - 1;
    for (int l = b; l < b + 2; l++)
    {
        // for each point, extra value for each derivative
        for (int i = 0; i <= (int)((n-1)/2); i++)
        {
            if (equidistant)
                wx.push_back(l);
            else
                wx.push_back(x[l]);
            wy.push_back(y[l]);
        }
    }

    //////std::cout << "wx" << std::endl;
    //////for (size_t i = 0; i < wx.size(); i++)
    //////    std::cout << wx[i] << " ";
    //////std::cout << std::endl;
    //////std::cout << "wy" << std::endl;
    //////for (size_t i = 0; i < wy.size(); i++)
    //////    std::cout << wy[i] << " ";
    //////std::cout << std::endl;


    // compute interpolation
    double retval = y[0];

    for (int j = 1; j <= n; j++)
    {
        std::vector<double> x1 = {};
        std::vector<double> y1 = {};
        for (int l = 0; l <= j; l++)
        {
            x1.push_back(wx[l]);
            y1.push_back(wy[l]);
            //std::cout << "Filling x1 with " << wx[l] << std::endl;
        }

        double d = diff_hermite(x1, y1, x_vec, y_vec, n, p, N);

        double prod = 1;
        for (int k = 0; k < j; k++)
            prod *= (x_sample - x1[k]);

        retval += d * prod;

        ///std::cout << "OUT for j = " << j << ": " << retval << std::endl << std::endl << std::endl;

    }

    ///std::cout << "END: " << retval << std::endl;

    return retval;
}



#endif
