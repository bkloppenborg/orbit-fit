/*
 * fitboth.cpp
 *
 *  Created on: Sep 6, 2011
 *      Author: bkloppenborg
 *
 *  Routines for fitting both RV and astrometric data simultaneously.
 *  Here we leverage existing routines in fitast_common and fitrv_common
 */

#include <string>
#include <vector>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <float.h>

#include "common.h"
#include "constants.h"
#include "fitast_common.h"
#include "fitrv_common.h"
#include "multinest_inf.h"

// Externals for the data
extern int n_rv_data;
extern int n_ast_data;
extern vector< vector<double> > rv_data;
extern vector< vector<double> > ast_data;

// Externals for some orbital parameters:
extern double Omega_min;
extern double Omega_max;
extern double inc_min;
extern double inc_max;
extern double omega_min;
extern double omega_max;

extern double prior_omega;
extern double prior_e;
extern double prior_T;
extern double prior_tau;

// Booleans for additional fitting parameters:
extern bool fit_turbulence;
extern bool fit_motion;
extern bool fit_astrometric_noise;

extern int motion_offset;

extern double fitrv::s_min;
extern double fitrv::s_max;
extern double fitast::s_min;
extern double fitast::s_max;

// Default number of parameters for the fitting routines.
int n_rv_params = 6;
int n_ast_params = 7;

using namespace std;

// Dumper
void dumper(int *nSamples, int *nlive, int *nPar, double **physLive, double **posterior, double *paramConstr, double *maxLogLike, double *logZ, double *logZerr)
{

}

// Log likelihood function, calls routines in fitrv and fitast
void log_likelihood(double *Cube, int *ndim, int *npars, double *lnew)
{
	// Because we use fitrv::log_likelihood and fitast::log_likelihood
	// we need to reroute the parameters.  We store them in these two local arrays
	// (local because of MPI)
	double * rv_params = new double[n_rv_params];
	double * ast_params = new double[n_ast_params];

	int ast_opt_offset = 8;

	// First reroute the parameters.  Indexing here is a little tricky because
	// some parameters are optional.  The default order is:
	// omega, e, tau, T, K, gamma, Omega, i, alpha, rv_s*, x0*, y0*, mu_x*, mu_y*, pi*, ast_s*
	// where the starred parameter are optional.
	rv_params[0] = Cube[0];	// omega
	rv_params[1] = Cube[1]; // e
	rv_params[2] = Cube[2]; // tau
	rv_params[3] = Cube[3]; // T
	rv_params[4] = Cube[4]; // K
	rv_params[5] = Cube[5]; // gamma

	if(fit_turbulence)
	{
		rv_params[6] = Cube[9]; // s (additional stellar jitter term)
		ast_opt_offset += 1;
	}

	ast_params[0] = Cube[0]; // omega
	ast_params[1] = Cube[1]; // e
	ast_params[2] = Cube[2]; // tau
	ast_params[3] = Cube[3]; // T
	ast_params[4] = Cube[6]; // Omega
	ast_params[5] = Cube[7]; // i
	ast_params[6] = Cube[8]; // alpha

	if(fit_motion)
	{
		ast_params[7]  = Cube[ast_opt_offset + 1]; // x0
		ast_params[8]  = Cube[ast_opt_offset + 2]; // y0
		ast_params[9]  = Cube[ast_opt_offset + 3]; // mu_x
		ast_params[10] = Cube[ast_opt_offset + 4]; // mu_y
		ast_params[11] = Cube[ast_opt_offset + 5]; // pi
	}

	if(fit_astrometric_noise)
	{
		ast_params[6 + motion_offset + 1] = Cube[ast_opt_offset + motion_offset + 1];
		ast_params[6 + motion_offset + 2] = Cube[ast_opt_offset + motion_offset + 2];
	}

	double lnew_ast = 0;
	double lnew_rv = 0;

	// Now call the fitting routines, have them compute the log likelihood.
	fitast::log_likelihood(ast_params, &n_ast_params, &n_ast_params, &lnew_ast);
	fitrv::log_likelihood(rv_params, &n_rv_params, &n_rv_params, &lnew_rv);

	// NOTE: Remember, we have double-counted the priors for omega, e, T, and tau so we need to
	// add that back in to not bias our result.  It's positive here.
	//Cube[0];	// omega
	//Cube[1]; // e
	//Cube[2]; // tau
	//Cube[3]; // T
	double d_cnt_priors = prior_omega
    		+ prior_e
    		+ 1.0 / Cube[3] * prior_T
    		+ 1.0 / Cube[4] * prior_tau;

	// Now, add the likelihoods together.
	*lnew = lnew_ast + lnew_rv + d_cnt_priors;

	//printf("%f %f %f\n", lnew_ast, lnew_rv, *lnew);

	// Now push the parameters back to multinest:
	// First reroute the parameters
	Cube[0] = rv_params[0];	// omega
	Cube[1] = rv_params[1]; // e
	Cube[2] = rv_params[2]; // tau
	Cube[3] = rv_params[3]; // T
	Cube[4] = rv_params[4]; // K
	Cube[5] = rv_params[5]; // gamma
	Cube[6] = ast_params[4]; // Omega
	Cube[7] = ast_params[5]; // i
	Cube[8] = ast_params[6]; // alpha

	if(fit_turbulence)
		Cube[9] = rv_params[6]; // s_rv (additional stellar jitter term)


	if(fit_motion)
	{
		Cube[ast_opt_offset + 1] = ast_params[7]; // x0
		Cube[ast_opt_offset + 2] = ast_params[8]; // y0
		Cube[ast_opt_offset + 3] = ast_params[9]; // mu_x
		Cube[ast_opt_offset + 4] = ast_params[10]; // mu_y
		Cube[ast_opt_offset + 5] = ast_params[11]; // pi
	}

	if(fit_astrometric_noise)
	{
		Cube[ast_opt_offset + motion_offset + 1] = ast_params[6 + motion_offset + 1];
		Cube[ast_opt_offset + motion_offset + 2] = ast_params[6 + motion_offset + 2];
	}

	// Free allocated memory:
	delete rv_params;
	delete ast_params;

}

// Prints out help describing the options on the command line
void print_help()
{
	string usage = "fitboth - A Bayesian Astrometric Orbit Fitting Program\n"
	"Fits Astrometric positions and Radial Velocities using MultiNest\n\n"
	"Usage: \n"
	" fitast rv_data_file astrometry_data_file [options] \n"
	"";

	cout << usage << "\n";
}

void run_fit()
{
	// Print out the limits on the parameters:
	printf("Starting fit with the following limits: \n");
	printf("Param       : Min        Max\n");
	print_common_param_limits();
	fitast::print_param_limits();
	fitrv::print_param_limits();

	fitast::compute_scales();
	fitrv::compute_scales();

	fitast::compute_partial_priors();
	fitrv::compute_partial_priors();

	n_rv_data = rv_data.size();
	n_ast_data = ast_data.size();

    printf("Found %i RV data points.\n", n_rv_data);
    printf("Found %i astrometry data points.\n", n_ast_data);

    // Determine the number of optional parameters and
    // allocate space for parameter sorting.
    if(fit_turbulence)
    	n_rv_params += 1;

    if(fit_motion)
    {
    	n_ast_params += 5;
    	motion_offset = 5;
    }

    if(fit_astrometric_noise)
    	n_ast_params += 2;


	// set the MultiNest sampling parameters
	int mmodal = 1;					// do mode separation?
	int ceff = 0;					// run in constant efficiency mode?
	int nlive = 1000;				// number of live points
	double efr = 1.0;				// set the required efficiency
	double tol = 0.5;				// tol, defines the stopping criteria
	int ndims = n_rv_params + n_ast_params - 4;					// dimensionality (no. of free parameters)
	int nPar = n_rv_params + n_ast_params - 4;					// total no. of parameters including free & derived parameters
	int nClsPar = n_rv_params + n_ast_params - 4;				// no. of parameters to do mode separation on
	int updInt = 100;				// after how many iterations feedback is required & the output files should be updated
									// note: posterior files are updated & dumper routine is called after every updInt*10 iterations
	double Ztol = -1E90;			// all the modes with logZ < Ztol are ignored
	int maxModes = 100;				// expected max no. of modes (used only for memory allocation)
	int pWrap[ndims];				// which parameters to have periodic boundary conditions?
	for(int i = 0; i < ndims; i++)
	    pWrap[i] = 0;

	// The angular parameters can have periodic boundary conditions:
	if(omega_min == 0 && omega_max == TWO_PI)
		pWrap[0] = 1;

	if(Omega_min == 0 && Omega_max == TWO_PI)
		pWrap[6] = 1;

	if(inc_min == -PI && inc_max == PI)
		pWrap[8] = 1;


	char root[100] = "chains/fitboth-";		// root for output files
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

// Parse command-line options that are specific to this program
void ParseProgOptions(int argc, char *argv[], bool & param_error)
{
	// Init values:

	for (int i = 1; i < argc; i++)
	{
		// Help
		if(strcmp(argv[i], "-h") == 0)
		{
			print_help();
			param_error = true;
		}

		if(strcmp(argv[i], "-s_min") == 0)
			param_error = true;

		if(strcmp(argv[i], "-s_max") == 0)
			param_error = true;


	}

	if(param_error)
		printf("Use of s_min or s_max is ambigious.  Use rv_s_min/max and ast_s_min/max intead.\n");
}

// The main routine.  Basically just used to parse out some parameters before handing
// things off to other functions.
int main(int argc, char *argv[])
{
	// Init/allocate locals:
	double tmp;
	bool param_error = false;

	// First parse command line arguments that are only for this program:
    if(argc == 1)
        print_help();

    if(argc < 3)
    	cout << "Missing filename on command line";


    // Read in the input filename:
    string input_rv = string(argv[1]);
    string input_ast = string(argv[2]);

    // Parse overrides
    ParseCommonParams(argc, argv, param_error);
    fitrv::ParseProgOptions(argc, argv, param_error);
    fitast::ParseProgOptions(argc, argv, param_error);
    ParseProgOptions(argc, argv, param_error);

    if(param_error)
    	return 0;

    // Read in the file of RV data.
    // Each row should contain time, RV, [errors]
    const string comment_chars("\\#~$&£%");

    vector< vector<int> > split_info;

    fitrv::read_data(input_rv, comment_chars, split_info);
    fitast::read_data(input_ast, comment_chars, split_info);

    run_fit();

	return 0;
}




