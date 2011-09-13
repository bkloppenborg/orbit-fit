/*
 * main.cpp
 *
 *  Created on: Aug 14, 2011
 *      Author: bkloppenborg
 *
 *  Fit radial velocity data using a nested sampling Bayesian routine.
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

#include "fitrv.h"
#include "ReadTextFile.h"
#include "orbit.h"
#include "common.h"
#include "constants.h"
#include "fitrv_common.h"

#include "multinest_inf.h"

using namespace std;
using namespace fitrv;

// Data and number of data.
namespace fitrv
{
	double default_error;
}
vector< vector<double> > rv_data;
int n_rv_data;

extern vector<string> param_names;

// Parameters for the fit:
double gamma_min;
double gamma_max;
double K_min;
double K_max;
extern double omega_min;
extern double omega_max;
extern double e_min;
extern double e_max;
extern double T_min;
extern double T_max;
extern double tau_min;
extern double tau_max;
namespace fitrv
{
	double s_min;
	double s_max;
}

// Scaling factors
double scale_gamma;
double scale_K;
extern double scale_omega;
extern double scale_e;
extern double scale_T;
extern double scale_tau;
namespace fitrv
{
	double scale_s;
}

// Partially (or fully) computed priors.
double prior_gamma;
double prior_K;
extern double prior_omega;
extern double prior_e;
extern double prior_T;
extern double prior_tau;
namespace fitrv
{
	double prior_s;
}

// Following the discussion in Feroz (2011) the semi-amplitude and jitter terms have value 1.
namespace fitrv
{
	double s_0 = 1;
}
double K_0 = 1;

// Counter for optional parameters:
namespace fitrv
{
	int opt_params = 0;
}
bool fit_turbulence = false;


void fitrv::read_data(string filename, string comment_chars, vector< vector<int> > split_info)
{
	// First determine the type of observation and fork it off to the
    vector < string > lines = ReadFile(filename, comment_chars, "Cannot Open Radial Velocity Data File");
	vector < vector<double> > results;

	int obs_type = -1;
	double t, rv, e_rv;

	for (unsigned int i = 0; i < lines.size(); i++)
	{
		vector < string > tokens;

		if(split_info.size() > 0)
			tokens = Tokenize(lines[i], split_info);
		else
			tokens = Tokenize(lines[i]);


		// And now attempt to read in the line
		try
		{
			t = atof(tokens[0].c_str());
			rv = atof(tokens[1].c_str());
		}
		catch(...)
		{
			throw std::runtime_error("Could not parse line in RV data file.");
		}

		// Now do the RV, permit
		try
		{
			e_rv = atof(tokens[2].c_str());
		}
		catch(...)
		{
			e_rv = 0;
		}

		if(e_rv == 0)
			e_rv = default_error;

		// Enable if you want to see the data.
		//printf("%f %f %f \n", t, rv, e_rv);

		// Push this station on to the list of stations for this array.
		vector<double> temp;
		temp.push_back(t);
		temp.push_back(rv);
		temp.push_back(e_rv);
		rv_data.push_back(temp);
	}
}

void fitrv::log_likelihood(double *Cube, int *ndim, int *npars, double *lnew)
{
	// Local variables
	double s = 0;
	double M, E;
	double t, rvi, e_rvi;
	double rv, err, tmp, cos_E, sin_E;

	double omega 	= Cube[0] * scale_omega + omega_min;
	double e 		= Cube[1] * scale_e + e_min;
	double tau 		= Cube[2] * scale_tau + tau_min;
	double T 		= Cube[3] * scale_T + T_min;
	double K 		= Cube[4] * scale_K + K_min;
	double gamma 	= Cube[5] * scale_gamma + gamma_min;

    if(fit_turbulence)
    	s = Cube[6] * scale_s + s_min;

    // Now set the parameters for feedback to multinest
    Cube[0] = omega * RAD_TO_DEG;
    Cube[1] = e;
    Cube[2] = tau;
    Cube[3] = T;
    Cube[4] = K;
    Cube[5] = gamma;

    if(fit_turbulence)
    	Cube[6] = s;

    // A few pre-computed values
    double beta = sqrt(1 - e*e);
    double N = ComputeN(T);
    double cos_omega = cos(omega);
    double sin_omega = sin(omega);

    // Note, priors are derived from scale factors (positive numbers) in most cases,
    // so here we need to multiply by -1 on prior values to get the right scaling.
    double prior = 1.0 / (K + K_0) * 1.0 / log(1 + (K_max / K_0)*pow(T_min/T, 1.0 / 3)*(1.0 / beta));
    		- prior_omega
    		- prior_e
    		- 1.0 / T * prior_T
    		- 1.0 / tau * prior_tau
    		- 1.0 / prior_gamma;

    if(fit_turbulence)
    	prior -= 1.0 / (s + s_0) * prior_s;

    double llike = -0.5 * n_rv_data * log(TWO_PI);

    for(register int i = 0; i < n_rv_data; i++)
    {
    	// Look up the data
    	t = rv_data[i][0];
    	rvi = rv_data[i][1];
    	e_rvi = rv_data[i][2];

    	tmp = e_rvi * e_rvi + s * s;

    	// The following is based on GetRV_K in the orbital_motion library.  We don't call the
    	// function to avoid a few small computational costs (like computing eta time and time again).
    	M = ComputeM(tau, N, t);
    	E = ComputeE(M, e);
        cos_E = cos(E);
        sin_E = sin(E);

    	rv = K * beta / (1 - e * cos_E) * (beta * cos_E * cos_omega - sin_E * sin_omega);
    	err = gamma + rv - rvi;

    	llike -= 0.5 * log(tmp) + err * err / (2.0 * tmp);
    }

    // printf("%0.4f %0.4f %0.4f %0.4f %0.4f %0.4f \n", K, omega * RAD_TO_DEG, e, T, tau, gamma);
    //printf("RV: %f %f\n", llike, prior);

	*lnew = llike + prior;
}

// Dumper
void fitrv::dumper(int *nSamples, int *nlive, int *nPar, double **physLive, double **posterior, double *paramConstr, double *maxLogLike, double *logZ, double *logZerr)
{

/*
//	 paramConstr(4*nPar):
//   paramConstr(1) to paramConstr(nPar)	     	= mean values of the parameters
//   paramConstr(nPar+1) to paramConstr(2*nPar)    	= standard deviation of the parameters
//   paramConstr(nPar*2+1) to paramConstr(3*nPar)  	= best-fit (maxlike) parameters
//   paramConstr(nPar*4+1) to paramConstr(4*nPar)  	= MAP (maximum-a-posteriori) parameters
*/
//	printf("npar %i\n", *nPar);
//	printf("Pointer %p\n", paramConstr);
//
//	printf("maxLogLike %f, logZ %f, logZerr %f\n", 	*maxLogLike, *logZ, *logZerr);
//
//	for(int i = 0; i < *nPar + 1; i++)
//		printf("%s: %e %e %e\n", param_names[i].c_str(), paramConstr[i], paramConstr[2* (*nPar) + i], paramConstr[3* (*nPar) + i]);
}


void fitrv::compute_scales()
{
	// All of the parameters have been set, compute scale factors:
	scale_K = K_max - K_min;
	scale_omega = omega_max - omega_min;
	scale_e = e_max - e_min;
	scale_tau = tau_max - tau_min;
	scale_T = T_max - T_min;
	scale_gamma = gamma_max - gamma_min;
	scale_s = s_max - s_min;
}

void fitrv::compute_partial_priors()
{
	// Now compute the (sometimes partial) priors:
	prior_K = 1.0 / scale_K;
	prior_omega = 1.0 / scale_omega;
	prior_e = 1.0 / scale_e;
	prior_tau = 1.0 / log(tau_max / tau_min);
	prior_T = 1.0 / log(T_max / T_min);
	prior_gamma = 1.0 / scale_gamma;
	prior_s = 1.0 / scale_s;
}

void fitrv::run_fit()
{
	// Print out what parameters are being used here:
	printf("Starting fit with the following limits: \n");
	printf("Param       : Min        Max\n");
	print_common_param_limits();
	fitrv::print_param_limits();

	compute_scales();
	compute_partial_priors();

    n_rv_data = rv_data.size();
    printf("Found %i RV data points.\n", n_rv_data);

    if(fit_turbulence)
    	opt_params += 1;


	// set the MultiNest sampling parameters
	int mmodal = 1;					// do mode separation?
	int ceff = 0;					// run in constant efficiency mode?
	int nlive = 1000;				// number of live points
	double efr = 1.0;				// set the required efficiency
	double tol = 0.5;				// tol, defines the stopping criteria
	int ndims = 6 + opt_params;					// dimensionality (no. of free parameters)
	int nPar = 6 + opt_params;					// total no. of parameters including free & derived parameters
	int nClsPar = 6 + opt_params;				// no. of parameters to do mode separation on
	int updInt = 100;				// after how many iterations feedback is required & the output files should be updated
									// note: posterior files are updated & dumper routine is called after every updInt*10 iterations
	double Ztol = -1E90;			// all the modes with logZ < Ztol are ignored
	int maxModes = 100;				// expected max no. of modes (used only for memory allocation)
	int pWrap[ndims];				// which parameters to have periodic boundary conditions?
	for(int i = 0; i < ndims; i++)
	    pWrap[i] = 0;

	// omega can have periodic boundary conditions if it occupies 0 ... 2pi
	if(omega_min == 0 && omega_max == TWO_PI)
		pWrap[0] = 1;

	char root[100] = "chains/fitrv-";		// root for output files
	int seed = -1;					// random no. generator seed, if < 0 then take the seed from system clock
	int fb = 1;					    // need feedback on standard output?
	int resume = 0;					// resume from a previous job?
	int outfile = 1;				// write output files?
	int initMPI = 1;				// initialize MPI routines?, relevant only if compiling with MPI
							        // set it to F if you want your main program to handle MPI initialization

	double logZero = -DBL_MAX;		// points with loglike < logZero will be ignored by MultiNest
	int context = 0;				// not required by MultiNest, any additional information user wants to pass



	// Call MultiNest:
	run_multinest(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero,
	log_likelihood, dumper, context);
}

// Parse parameters not handled in ParseCommonParams
void fitrv::ParseProgOptions(int argc, char *argv[], bool & param_error)
{
	K_min = 0;
	K_max = 100;
	gamma_min = 10;
	gamma_max = -10;
	s_min = 0;
	s_max = 10;

	for (int i = 1; i < argc; i++)
	{
		// First see if the user is requesting help:
		if(strcmp(argv[i], "-h") == 0)
		{
			print_help();
			param_error = true;
		}

		// Fit atmospheric turbulence, treated as gaussian noise.
		if(strcmp(argv[i], "-turb") == 0)
		{
			printf("NOTE: Including atmospheric turbulence in fit.\n");
			fit_turbulence = true;
		}

		if(strcmp(argv[i], "-K_min") == 0)
			K_min = atof(argv[i + 1]);

		if(strcmp(argv[i], "-K_max") == 0)
			K_max = atof(argv[i + 1]);

		if(strcmp(argv[i], "-gamma_min") == 0)
			gamma_min = atof(argv[i + 1]);

		if(strcmp(argv[i], "-gamma_max") == 0)
			gamma_max = atof(argv[i + 1]);

		if(strcmp(argv[i], "-s_min") == 0 || strcmp(argv[i], "-rv_s_min") == 0)
			s_min = atof(argv[i + 1]);

		if(strcmp(argv[i], "-s_max") == 0 || strcmp(argv[i], "-rv_s_max") == 0)
			s_max = atof(argv[i + 1]);

		if(strcmp(argv[i], "-rv_err") == 0)
			default_error = atof(argv[i + 1]);
	}

	if(gamma_min == 0)
		gamma_min = -K_max;
	if(gamma_max == 0)
		gamma_max = K_max;
}

void fitrv::print_param_limits()
{
	printf("K (km/s)    : %1.4e %1.4e \n", K_min, K_max);
	printf("gamma (km/s): %1.4e %1.4e \n", gamma_min, gamma_max);

	if(fit_turbulence)
		printf("s (km/s)    : %1.4e %1.4e \n", s_min, s_max);
}
