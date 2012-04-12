/*
 * fitast_common.cpp
 *
 *  Created on: Aug 14, 2011
 *      Author: bkloppenborg
 *
 *  Fit astrometric data using a nested sampling Bayesian algorithm.
 */

#include <string>
#include <cstring>
#include <iostream>
#include <cstdio>
#include <stdexcept>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <float.h>

#include "fitast.h"
#include "ReadTextFile.h"
#include "orbit.h"
#include "common.h"
#include "constants.h"

#include "multinest.h"

using namespace std;
using namespace fitast;

// Global Variables (scales and (sometimes partial) priors:
double x_0_min;
double x_0_max;
double y_0_min;
double y_0_max;
double mu_x_min;
double mu_x_max;
double mu_y_min;
double mu_y_max;
double pi_min;
double pi_max;
namespace fitast
{
	double s_min;
	double s_max;
}

double Omega_min;
double Omega_max;
double inc_min;
double inc_max;
extern double omega_min;
extern double omega_max;
double alpha_min;
double alpha_max;
extern double e_min;
extern double e_max;
extern double tau_min;
extern double tau_max;
extern double T_min;
extern double T_max;

double scale_x_0;
double scale_y_0;
double scale_mu_x;
double scale_mu_y;
double scale_pi;
namespace fitast
{
	double scale_s;
}

double scale_Omega;
double scale_inc;
extern double scale_omega;
double scale_alpha;
extern double scale_e;
extern double scale_tau;
extern double scale_T;

double prior_x_0;
double prior_y_0;
double prior_mu_x;
double prior_mu_y;
double prior_pi;
double prior_Omega;
double prior_inc;
extern double prior_omega;
double prior_alpha;
extern double prior_e;
extern double prior_tau;
extern double prior_T;
namespace fitast
{
	double prior_s;
}

extern string output;

vector< vector<double> > ast_data;
int n_ast_data;

// Optional parameters
// orbit_param_offset = 0 -> don't fit zero point, proper motion, parallax, etc.  Other valid value = 5.
namespace fitast
{
	int opt_params = 0;
}

int motion_offset = 0;

bool fit_motion = false;

namespace fitast
{
	bool read_no_error = false;
	double default_error;
}
bool fit_astrometric_noise = false;


void fitast::log_likelihood(double * params, int & ndim, int & npars, double & lnew, void * misc)
{
	// Locals
	double t, xi, e_xi, yi, e_yi, P_a, P_d;
	double x, y, err_x, err_y;
	double M, E, cos_E, sin_E;
	double x_0, y_0, mu_x, mu_y, pi;
	double dt;
	double s_x, s_y;

	s_x = 0;
	s_y = 0;

	// Pull out the parameters from the cube
	double omega 	= params[0] * scale_omega + omega_min;
	double e 		= params[1] * scale_e + e_min;
	double tau 		= params[2] * scale_tau + tau_min;
	double T 		= params[3] * scale_T + T_min;
	double Omega 	= params[4] * scale_Omega + Omega_min;
	double inc 		= params[5] * scale_inc + inc_min;
	double alpha 	= params[6] * scale_alpha + alpha_min;

	if(fit_motion)
	{
		x_0 	= params[7] * scale_x_0 + x_0_min;
		y_0 	= params[8] * scale_y_0 + y_0_min;
		mu_x 	= params[9] * scale_mu_x + mu_x_min;
		mu_y 	= params[10] * scale_mu_y + mu_y_min;
		pi 		= params[11] * scale_pi + pi_min;
	}

	if(fit_astrometric_noise)
	{
		s_x = params[6 + motion_offset + 1] * scale_s + s_min;
		s_y = params[6 + motion_offset + 2] * scale_s + s_min;
	}

	// Now set the scaled parameters back in the cube:
	params[0] = omega * RAD_TO_DEG;
	params[1] = e;
	params[2] = tau;
	params[3] = T;
	params[4] = Omega * RAD_TO_DEG;
	params[5] = inc * RAD_TO_DEG;
	params[6] = alpha;

	if(fit_motion)
	{
		params[7] = x_0;
		params[8] = y_0;
		params[9] = mu_x;
		params[10] = mu_y;
		params[11] = pi;
	}

	if(fit_astrometric_noise)
	{
		params[6 + motion_offset + 1] = s_x;
		params[6 + motion_offset + 2] = s_y;
		//printf("s_x %i %e s_y %i %e\n", 6 + motion_offset + 1, s_x, 6 + motion_offset + 2, s_y);
	}

    // Note, if an error is found here, be sure to update the GetRV and GetXY functions.
    double m1 = 0;
    double l1 = 0;
    double m2 = 0;
    double l2 = 0;
    double junk = 0;
    Compute_Coefficients(Omega, inc, omega, l1, m1, junk, l2, m2, junk);
    
    // A few pre-computed values
    double beta = sqrt(1 - e*e);
    double n = ComputeN(T);
    double s_x2 = s_x * s_x;
    double s_y2 = s_y * s_y;

    // Compute the prior:
    double prior = prior_omega
    		+ prior_e
    		+ 1.0 / tau * prior_tau
    		+ 1.0 / T * prior_T
    		+ prior_Omega
    		+ prior_inc
    		+ prior_alpha;

    if(fit_motion)
    	prior += prior_x_0 + prior_y_0 + prior_mu_x + prior_mu_y + prior_pi;

    if(fit_astrometric_noise)
    	prior += 2*prior_s;

    double llike = -n_ast_data * log(TWO_PI);

    for(register int i = 0; i < n_ast_data; i++)
    {
    	// Pull out the data:
    	t = ast_data[i][0];
    	xi = ast_data[i][1];
    	e_xi = ast_data[i][2];
    	yi = ast_data[i][3];
    	e_yi = ast_data[i][4];
    	P_a = ast_data[i][5];
    	P_d = ast_data[i][6];

    	// Compute orbital elements:
        M = ComputeM(tau, n, t);
        E = ComputeE(M, e);
        cos_E = cos(E);
        sin_E = sin(E);

    	// Get the XY positions of the orbit
    	Compute_xy(alpha, beta, e, l1, l2, m1, m2, cos_E, sin_E, x, y);

    	if(fit_motion)
    	{
    		dt = t - tau;
    		x += x_0 + mu_x * dt + pi * P_a;
    		y += y_0 + mu_y * dt + pi * P_d;
    	}

    	err_x = x - xi;
    	err_y = y - yi;

    	llike -= log((e_xi + s_x) * (e_yi + s_y))
    			+ 0.5 * (
    				  err_x * err_x / (e_xi * e_xi + s_x2)
    				+ err_y * err_y / (e_yi * e_yi + s_y2)
    			);

    }

	//printf("%e %f %e %f\n", mu_x, dt, mu_x * dt, 0.0);
    //printf("%f %f %f %f %f %f \n", x, xi, e_xi, y, yi, e_yi);
    //printf("%e %e %e %e\n", x, xi, err_x, err_x * err_x / (2 * (e_xi * e_xi + s_x2)));
    //printf("AST: %f %f\n", llike, prior);

    // Assign the value and we're done.
	lnew = llike + prior;
}

void fitast::dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double ** paramConstr, double &maxLogLike, double & logZ, double & logZerr, void * misc)
{
	// do nothing for now.
}


void fitast::compute_scales()
{
	// All of the parameters have been set, compute scale factors:
	scale_Omega = Omega_max - Omega_min;
	scale_inc = inc_max - inc_min;
	scale_omega = omega_max - omega_min;
	scale_alpha = alpha_max - alpha_min;
	scale_e = e_max - e_min;
	scale_tau = tau_max - tau_min;
	scale_T = T_max - T_min;
	scale_x_0 = x_0_max - x_0_min;
	scale_y_0 = y_0_max - y_0_min;
	scale_mu_x = mu_x_max - mu_x_min;
	scale_mu_y = mu_y_max - mu_x_min;
	scale_pi = pi_max - pi_min;
	scale_s = s_max - s_min;
}

void fitast::compute_partial_priors()
{
	// Now compute the (sometimes partial) priors:
	// NOTE: Some of the priors are negative because the scales are max-min instead of min-max.
	prior_Omega = -1.0 / scale_Omega;
	prior_inc = -1.0 / scale_inc;
	prior_omega = -1.0 / scale_omega;
	prior_alpha = -1.0 / scale_alpha;
	prior_e = -1.0 / scale_e;
	prior_tau = 1.0 / log(tau_max / tau_min);
	prior_T = 1.0 / log(T_max / T_min);
	prior_x_0 = -1.0 / scale_x_0;
	prior_y_0 = -1.0 / scale_y_0;
	prior_mu_x = -1.0 / scale_mu_x;
	prior_mu_y = -1.0 / scale_mu_y;
	prior_pi = -1.0 / scale_pi;
	prior_s = -1.0 / scale_s;
}

void fitast::run_fit()
{
	// Print out what parameters are being used here:
	printf("Starting fit with the following limits: \n");
	printf("Param       : Min        Max\n");
	print_common_param_limits();
	fitast::print_param_limits();

	compute_scales();
	compute_partial_priors();

    n_ast_data = ast_data.size();
    printf("Found %i data astrometric points.\n", n_ast_data);
    opt_params = 0;

    if(fit_motion)
    {
    	opt_params += 5;
    	motion_offset = 5;
    }

    if(fit_astrometric_noise)
    	opt_params += 2;

	// set the MultiNest sampling parameters
	int mmodal = 1;					// do mode separation?
	int ceff = 0;					// run in constant efficiency mode?
	int nlive = 1000;				// number of live points
	double efr = 1.0;				// set the required efficiency
	double tol = 0.5;				// tol, defines the stopping criteria
	int ndims = 7 + opt_params;		// dimensionality (no. of free parameters)
	int nPar = 7 + opt_params;		// total no. of parameters including free & derived parameters
	int nClsPar = 7 + opt_params;	// no. of parameters to do mode separation on
	int updInt = 100;				// after how many iterations feedback is required & the output files should be updated
									// note: posterior files are updated & dumper routine is called after every updInt*10 iterations
	double Ztol = -1E90;			// all the modes with logZ < Ztol are ignored
	int maxModes = 100;				// expected max no. of modes (used only for memory allocation)
	int pWrap[ndims];				// which parameters to have periodic boundary conditions?
	for(int i = 0; i < ndims; i++)
	    pWrap[i] = 0;

	// Enable wrapping for some parameters if they occupy the full range:
	if(omega_min == 0 && omega_max == TWO_PI)
		pWrap[0] = 1;

	if(Omega_min == 0 && Omega_max == TWO_PI)
		pWrap[4] = 1;

	if(inc_min == -PI && inc_max == PI)
		pWrap[5] = 1;


	int seed = -1;					// random no. generator seed, if < 0 then take the seed from system clock
	int fb = 1;					    // need feedback on standard output?
	int resume = 1;					// resume from a previous job?
	int outfile = 1;				// write output files?
	int initMPI = 0;				// initialize MPI routines?, relevant only if compiling with MPI
							        // set it to F if you want your main program to handle MPI initialization

	double logZero = -DBL_MAX;		// points with loglike < logZero will be ignored by MultiNest
//	int context = 0;				// not required by MultiNest, any additional information user wants to pass
	int maxIterations = 1E9;

	void * misc = NULL;

	string temp;
	if(output.length() > 0)
		temp = output;
	else
		temp = "chains/fitast-";

	const std::string path = temp;

    // Run the nested sampling algorithm
    nested::run(mmodal, ceff, nlive, tol,
        efr, ndims, nPar, nClsPar,
        maxModes, updInt, Ztol, path,
        seed, pWrap, fb, resume,
        outfile, initMPI, logZero, maxIterations,
        fitast::log_likelihood,
        fitast::dumper,
        misc);
}

// Parse command-line options that are specific to this program
void fitast::ParseProgOptions(int argc, char *argv[], bool & param_error)
{
	// Init values:
	default_error = 1;
	alpha_min = 0;
	alpha_max = 1;
	Omega_min = 0;
	Omega_max = TWO_PI;
	inc_min = -PI;
	inc_max = PI;

	x_0_min = 0;
	x_0_max = 180;
	y_0_min = -90;
	y_0_max = 90;
	mu_x_min = 0;
	mu_x_max = 1;
	mu_y_min = 0;
	mu_y_max = 1;
	pi_min = 0;
	pi_max = 1;
	s_min = 0;
	s_max = 1;

	for (int i = 1; i < argc; i++)
	{
		// Help
		if(strcmp(argv[i], "-h") == 0)
		{
			print_help();
			param_error = true;
		}

		// First see if the user is requesting help:
		if(strcmp(argv[i], "-motion") == 0)
			fit_motion = true;

//		if(strcmp(argv[i], "-r_theta") == 0)
//			read_r_theta = true;

		if(strcmp(argv[i], "-noerr") == 0)
			read_no_error = true;

		if(strcmp(argv[i], "-ast_noise") == 0)
			fit_astrometric_noise = true;

		if(strcmp(argv[i], "-mu_x_min") == 0)
			mu_x_min = atof(argv[i+1]);

		if(strcmp(argv[i], "-mu_x_max") == 0)
			mu_x_max = atof(argv[i+1]);

		if(strcmp(argv[i], "-mu_y_min") == 0)
			mu_y_min = atof(argv[i+1]);

		if(strcmp(argv[i], "-mu_y_max") == 0)
			mu_y_max = atof(argv[i+1]);

		if(strcmp(argv[i], "-pi_min") == 0)
			pi_min = atof(argv[i+1]);

		if(strcmp(argv[i], "-pi_max") == 0)
			pi_max = atof(argv[i+1]);

		if(strcmp(argv[i], "-x0_min") == 0)
			x_0_min = atof(argv[i+1]);

		if(strcmp(argv[i], "-x0_max") == 0)
			x_0_max = atof(argv[i+1]);

		if(strcmp(argv[i], "-y0_min") == 0)
			y_0_min = atof(argv[i+1]);

		if(strcmp(argv[i], "-y0_max") == 0)
			y_0_max = atof(argv[i+1]);

		if(strcmp(argv[i], "-Omega_min") == 0)
			Omega_min = atof(argv[i + 1]) * DEG_TO_RAD;

		if(strcmp(argv[i], "-Omega_max") == 0)
			Omega_max = atof(argv[i + 1]) * DEG_TO_RAD;

		if(strcmp(argv[i], "-inc_min") == 0)
			inc_min = atof(argv[i + 1]) * DEG_TO_RAD;

		if(strcmp(argv[i], "-inc_max") == 0)
			inc_max = atof(argv[i + 1]) * DEG_TO_RAD;

		if(strcmp(argv[i], "-alpha_min") == 0)
			alpha_min = atof(argv[i + 1]);

		if(strcmp(argv[i], "-alpha_max") == 0)
			alpha_max = atof(argv[i + 1]);

		if(strcmp(argv[i], "-ast_err") == 0)
			default_error = atof(argv[i + 1]);

		if(strcmp(argv[i], "-s_min") == 0 || strcmp(argv[i], "-ast_s_min") == 0)
			s_min = atof(argv[i + 1]);

		if(strcmp(argv[i], "-s_max") == 0 || strcmp(argv[i], "-ast_s_max") == 0)
			s_max = atof(argv[i + 1]);

    }

	if(fit_motion)
		printf("NOTE: Fitting zero points, proper motions, and parallax.\n");

//	if(read_r_theta)
//		printf("NOTE: Reading in data in r-theta format.  Converting to (x,y).\n");

	if(read_no_error)
		printf("NOTE: Reading data file without error columns.  \n"
			   "      Specify -ast_err n.nn to indicate the default error otherwise 1.00 is used.\n");
	if(default_error != 1)
		printf("NOTE: Using user-specified astrometric error of %e \n", default_error);

	if(fit_astrometric_noise)
		printf("NOTE: Including additional noise source (Gaussian) in astrometric fitting.\n");
}

void fitast::print_param_limits()
{
	// Setup the interface to multinest, run it.
	// Print out what parameters are being used here:
	printf("Omega (deg) : %1.4e %1.4e \n", Omega_min * RAD_TO_DEG, Omega_max * RAD_TO_DEG);
	printf("inc   (deg) : %1.4e %1.4e \n", inc_min * RAD_TO_DEG, inc_max * RAD_TO_DEG);
	printf("alpha (arb) : %1.4e %1.4e \n", alpha_min, alpha_max);

	if(fit_motion)
	{
		printf("x0 (arb u)  : %1.4e %1.4e \n", x_0_min, x_0_max);
		printf("y0 (arb u)  : %1.4e %1.4e \n", y_0_min, y_0_max);
		printf("mu_x (u/t)  : %1.4e %1.4e \n", mu_x_min, mu_x_max);
		printf("mu_y (u/t)  : %1.4e %1.4e \n", mu_y_min, mu_y_max);
		printf("pi (arb u)  : %1.4e %1.4e \n", pi_min, pi_max);
	}

	if(fit_astrometric_noise)
	{
		printf("s_(x,y)     : %1.4e %1.4e \n", s_min, s_max);
	}
}
