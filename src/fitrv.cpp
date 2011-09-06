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

#include "multinest_inf.h"

using namespace std;

// Data and number of data.
vector< vector<double> > data;
int n_data;

vector<string> param_names;

// Parameters for the fit:
double default_error;
double gamma_min;
double gamma_max;
double K_min;
double K_max;
double e_min;
double e_max;
double omega_min;
double omega_max;
double T_min;
double T_max;
double tau_min;
double tau_max;
double s_min;
double s_max;

// Scaling factors
double scale_gamma;
double scale_K;
double scale_e;
double scale_omega;
double scale_T;
double scale_tau;
double scale_s;

// Partially (or fully) computed priors.
double prior_gamma;
double prior_K;
double prior_e;
double prior_omega;
double prior_T;
double prior_tau;
double prior_s;

// Following the discussion in Feroz (2011) the semi-amplitude and jitter terms have value 1.
double s_0 = 1;
double K_0 = 1;

// Counter for optional parameters:
int opt_params = 0;
bool fit_turbulence = false;

// Prints out help describing the options on the command line
void print_help()
{
	string usage = "fitrv - A Bayesian Radial Velocity Fitting Program\n"
	"Fits Radial Velocity data using MultiNest\n\n"
	"Usage: \n"
	" fitrv input_file \n"
	" \n"
	"Optional Arguments: \n"
	" -h        Prints this message \n"
	" -turb     Include Atmospheric Turbulence component as a Gaussian \n"
	"           in the fitting routines. \n"
	"Overriding Limits: \n"
	" gamma, K, omega, e, tau, T, s \n"
	"  -*_min    Override lower bound on above parameter \n"
	"  -*_max    Override upper bound on above parameter \n"
	"\n"
	"   gamma, K, s: (km/s) \n"
	"   omega (deg) \n"
	"   e (unitless, [0, 1]\n"
	"   tau, T (time, arbitrary)";

	cout << usage << "\n";
}

void read_data(string filename, string comment_chars, double default_error, vector< vector<int> > split_info, vector< vector<double> > & data)
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
		data.push_back(temp);
	}
}

void log_likelihood(double *Cube, int *ndim, int *npars, double *lnew)
{
	// Local variables
	double s = 0;
	double M, E;
	double t, rvi, e_rvi;
	double rv, err, tmp, cos_E, sin_E;

	double K = Cube[0] * scale_K + K_min;
	double omega = Cube[1] * scale_omega + omega_min;
	double e = Cube[2] * scale_e + e_min;
	double tau = Cube[3] * scale_tau + tau_min;
	double T = Cube[4] * scale_T + T_min;
	double gamma = Cube[5] * scale_gamma + gamma_min;

    if(fit_turbulence)
    	s = Cube[6] * scale_s + s_min;

    // Now set the parameters for feedback to multinest
    Cube[0] = K;
    Cube[1] = omega * RAD_TO_DEG;
    Cube[2] = e;
    Cube[3] = tau;
    Cube[4] = T;
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
    double prior = 1.0 / (K + K_0) * 1.0 / log(1 + (K_max / K_0)*pow(T_min/T, 1.0 / 3)*(1.0 / beta))
    		- prior_omega
    		- prior_e
    		- 1.0 / T * prior_T
    		- 1.0 / tau * prior_tau
    		- 1.0 / prior_gamma;

    if(fit_turbulence)
    	prior -= 1.0 / (s + s_0) * prior_s;

    double llike = -0.5 * n_data * log(TWO_PI);

    for(register int i = 0; i < n_data; i++)
    {
    	// Look up the data
    	t = data[i][0];
    	rvi = data[i][1];
    	e_rvi = data[i][2];

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
    // printf("%0.4f %0.4f %0.4f\n", rv, rvi, llike);

	*lnew = llike + prior;
}

// Dumper
void dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double *paramConstr, double &maxLogLike, double &logZ, double &logZerr)
{

/*
//	 paramConstr(4*nPar):
//   paramConstr(1) to paramConstr(nPar)	     	= mean values of the parameters
//   paramConstr(nPar+1) to paramConstr(2*nPar)    	= standard deviation of the parameters
//   paramConstr(nPar*2+1) to paramConstr(3*nPar)  	= best-fit (maxlike) parameters
//   paramConstr(nPar*4+1) to paramConstr(4*nPar)  	= MAP (maximum-a-posteriori) parameters
*/
//	printf("npar %i\n", nPar);
//
//	for(register int i = 0; i < nPar; i++)
//		printf("%s %1.4e %1.4e %1.4e %1.4e\n", param_names[i].c_str(), paramConstr[i], 0.0, 0.0, 0.0);

}



void run_fit(vector< vector<double> > & data)
{
	// Print out what parameters are being used here:
	printf("Starting fit with the following limits: \n");
	printf("Param       : Min        Max\n");
	printf("K (km/s)    : %1.4e %1.4e \n", K_min, K_max);
	printf("omega (deg) : %1.4e %1.4e \n", omega_min * RAD_TO_DEG, omega_max * RAD_TO_DEG);
	printf("e           : %1.4e %1.4e \n", e_min, e_max);
	printf("T (time)    : %1.4e %1.4e \n", T_min, T_max);
	printf("tau (time)  : %1.4e %1.4e \n", tau_min, tau_max);
	printf("gamma (km/s): %1.4e %1.4e \n", gamma_min, gamma_max);

	if(fit_turbulence)
		printf("s (km/s)    : %1.4e %1.4e \n", s_min, s_max);

	// All of the parameters have been set, compute scale factors:
	scale_K = K_max - K_min;
	scale_omega = omega_max - omega_min;
	scale_e = e_max - e_min;
	scale_tau = tau_max - tau_min;
	scale_T = T_max - T_min;
	scale_gamma = gamma_max - gamma_min;
	scale_s = s_max - s_min;

	// Now compute the (sometimes partial) priors:
	prior_K = 1.0 / scale_K;
	prior_omega = 1.0 / scale_omega;
	prior_e = 1.0 / scale_e;
	prior_tau = 1.0 / log(tau_max / tau_min);
	prior_T = 1.0 / log(T_max / T_min);
	prior_gamma = 1.0 / scale_gamma;
	prior_s = 1.0 / scale_s;

    n_data = data.size();
    printf("Found %i data points.\n", n_data);

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
		pWrap[1] = 1;

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
void ParseProgOptions(int argc, char *argv[], bool & param_error)
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

		if(strcmp(argv[i], "-s_min") == 0)
			s_min = atof(argv[i + 1]);

		if(strcmp(argv[i], "-s_max") == 0)
			s_max = atof(argv[i + 1]);

		if(strcmp(argv[i], "-err") == 0)
			default_error = atof(argv[i + 1]);
	}

	if(gamma_min == 0)
		gamma_min = -K_max;
	if(gamma_max == 0)
		gamma_max = K_max;
}

// The main routine.  Basically just used to parse out some parameters before handing
// things off to other functions.
int main(int argc, char *argv[])
{
	double tmp;
	bool param_error = false;

    if(argc == 1)
    {
        print_help();
        return 0;
    }

    if(argc < 2)
    	cout << "Missing filename on command line.";

    // Read in the input filename:
    string input_rv = string(argv[1]);

    // Parse the remaining command-line options
    ParseCommonParams(argc, argv, omega_min, omega_max, e_min, e_max, tau_min, tau_max, T_min, T_max, param_error);
    ParseProgOptions(argc, argv, param_error);

    // If there was an error, quit.
    if(param_error)
    	return 0;

    // Read in the file of RV data.
    // Each row should contain time, RV, [errors]
    const string comment_chars("\\#~$&£%");

    int pos_size[] = {0, 12, 13, 6, 20, 4};
    vector< vector<int> > split_info;
//    for(int i = 0; i < 3; i++)
//    {
//    	vector<int> tmp;
//    	tmp.push_back(pos_size[2*i]);
//    	tmp.push_back(pos_size[2*i+1]);
//    	split_info.push_back(tmp);
//    }

    // Push the parameter names onto the vector.
    param_names.push_back("K     ");
    param_names.push_back("omega ");
    param_names.push_back("e     ");
    param_names.push_back("T     ");
    param_names.push_back("tau   ");
    param_names.push_back("gamma ");
    param_names.push_back("s     ");

    read_data(input_rv, comment_chars, default_error, split_info, data);

    run_fit(data);

	return 0;
}
