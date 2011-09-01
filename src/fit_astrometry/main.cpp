/*
 * main.cpp
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

#include "main.h"
#include "ReadTextFile.h"
#include "orbit.h"
#include "common.h"

#include "multinest_inf.h"

using namespace std;

// Global Variables (scales and (sometimes partial) priors:
double default_error;
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

double Omega_min;
double Omega_max;
double inc_min;
double inc_max;
double omega_min;
double omega_max;
double alpha_min;
double alpha_max;
double e_min;
double e_max;
double tau_min;
double tau_max;
double T_min;
double T_max;

double scale_x_0;
double scale_y_0;
double scale_mu_x;
double scale_mu_y;
double scale_pi;

double scale_Omega;
double scale_inc;
double scale_omega;
double scale_alpha;
double scale_e;
double scale_tau;
double scale_T;

double prior_x_0;
double prior_y_0;
double prior_mu_x;
double prior_mu_y;
double prior_pi;
double prior_Omega;
double prior_inc;
double prior_omega;
double prior_alpha;
double prior_e;
double prior_tau;
double prior_T;

vector< vector<double> > data;
int n_data;

// Optional parameters
// orbit_param_offset = 0 -> don't fit zero point, proper motion, parallax, etc.  Other valid value = 5.
int opt_params = 0;

bool fit_motion = false;


// Prints out help describing the options on the command line
void print_help()
{
	string usage = "fitast - A Bayesian Astrometric Orbit Fitting Program\n"
	"Fits Astrometric Position data using MultiNest\n\n"
	"Usage: \n"
	" fitast input_file -units [rad|deg|mas] ... \n"
	" \n"
	"Optional Arguments: \n"
	" -h        Prints this message \n"
	" -motion   Fits zero points, proper motion, parallax in addition\n"
	"           to normal orbital parameters.\n"
	"Overriding Limits: \n"
	" Omega, inc, omega, alpha, e, tau, T \n"
	" -*_min    Override lower bound on above parameter \n"
	" -*_max    Override upper bound on above parameter \n"
	"";

	cout << usage << "\n";
}

void read_data(string filename, string comment_chars, double default_error, vector< vector<int> > split_info, vector< vector<double> > & data)
{
	// First determine the type of observation and fork it off to the
    vector < string > lines = ReadFile(filename, comment_chars, "Cannot Open Astrometry Data File");
	vector < vector<double> > results;

	double t, x, y, e_x, e_y, P_alpha, P_delta;

	// Now iterate through the lines and tokenize them.  Notice, we use vector.at(n) instead of vector[n]
	// so that we can catch the signals from exceptions.  The data file is only parsed once so this is ok.

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
			t = atof(tokens.at(0).c_str());
			x = atof(tokens.at(1).c_str());
			y = atof(tokens.at(3).c_str());
		}
		catch(...)
		{
			throw std::runtime_error("Could not parse line in astrometric data file.");
		}

		// Now for the uncertainties (these may not exist)
		try
		{
			e_x = atof(tokens.at(2).c_str());
		}
		catch(...)
		{
			e_x = 0;
		}

		try
		{
			e_y = atof(tokens.at(4).c_str());
		}
		catch(...)
		{
			e_y = 0;
		}

		// Lastly the parallax factors
		try
		{
			P_alpha = atof(tokens.at(5).c_str());
		}
		catch(...)
		{
			P_alpha = 0;
		}

		try
		{
			P_delta = atof(tokens.at(6).c_str());
		}
		catch(...)
		{
			P_delta = 0;
		}

		if(e_x == 0)
			e_x = default_error;
		if(e_y == 0)
			e_y = default_error;

		// Enable if you want to see the data.
		//printf("%f %f %f %f %f \n", t, x, e_x, y, e_y);

		// Push this station on to the list of stations for this array.
		vector<double> temp;
		temp.push_back(t);
		temp.push_back(x);
		temp.push_back(e_x);
		temp.push_back(y);
		temp.push_back(e_y);
		temp.push_back(P_alpha);
		temp.push_back(P_delta);
		data.push_back(temp);
	}
}

void log_likelihood(double *Cube, int *ndim, int *npars, double *lnew)
{
	// Locals
	double t, xi, e_xi, yi, e_yi, P_a, P_d;
	double x, y, err_x, err_y;
	double M, E, cos_E, sin_E;
	double x_0, y_0, mu_x, mu_y, pi;
	double dt;

	// Pull out the parameters from the cube
	double Omega = Cube[0] * scale_Omega + Omega_min;
	double inc = Cube[1] * scale_inc + inc_min;
	double omega = Cube[2] * scale_omega + omega_min;
	double alpha = Cube[3] * scale_alpha + alpha_min;
	double e = Cube[4] * scale_e + e_min;
	double tau = Cube[5] * scale_tau + tau_min;
	double T = Cube[6] * scale_T + T_min;

	if(fit_motion)
	{
		x_0 = Cube[7] * scale_x_0 + x_0_min;
		y_0 = Cube[8] * scale_y_0 + y_0_min;
		mu_x = Cube[9] * scale_mu_x + mu_x_min;
		mu_y = Cube[10] * scale_mu_y + mu_y_min;
		pi = Cube[11] * scale_pi + pi_min;
	}

	// Now set the scaled parameters back in the cube:
	Cube[0] = Omega * RAD_TO_DEG;
	Cube[1] = inc * RAD_TO_DEG;
    Cube[2] = omega * RAD_TO_DEG;
    Cube[3] = alpha;
    Cube[4] = e;
    Cube[5] = tau;
    Cube[6] = T;

	if(fit_motion)
	{
		Cube[7] = x_0;
		Cube[8] = y_0;
		Cube[9] = mu_x;
		Cube[10] = mu_y;
		Cube[11] = pi;
	}

    // Pre-compute a few values that are used frequently:
    double c_Omega = cos(Omega);
    double s_Omega = sin(Omega);
    double c_inc = cos(inc);
    double s_inc = sin(inc);
    double c_omega = cos(omega);
    double s_omega = sin(omega);

    // Note, if an error is found here, be sure to update the GetRV and GetXY functions.
    double l1 = c_Omega * c_omega - s_Omega * s_omega * c_inc;
    double m1 = s_Omega * c_omega + c_Omega * s_omega * c_inc;
    double l2 = -c_Omega * s_omega - s_Omega * c_omega * c_inc;
    double m2 = -s_Omega * s_omega + c_Omega * c_omega * c_inc;

    // A few pre-computed values
    double beta = sqrt(1 - e*e);
    double n = ComputeN(T);

    double prior = 0;

    if(fit_motion)
    	prior -= prior_x_0 + prior_y_0 + prior_mu_x + prior_mu_y + prior_pi;

    double llike = 0;

    for(register int i = 0; i < n_data; i++)
    {
    	// Pull out the data:
    	t = data[i][0];
    	xi = data[i][1];
    	e_xi = data[i][2];
    	yi = data[i][3];
    	e_yi = data[i][4];
    	P_a = data[i][5];
    	P_d = data[i][6];

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

    	llike -= log(TWO_PI * e_xi * e_yi)
    			+ err_x * err_x / (2 * e_xi * e_xi)
    			+ err_y * err_y / (2 * e_yi * e_yi);
    }

    //printf("%f %f %f %f %f %f \n", x, xi, e_xi, y, yi, e_yi);

	// Assign the value and we're done.
	*lnew = llike + prior;
}

void dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double *paramConstr, double &maxLogLike, double &logZ, double &logZerr)
{
	// do nothing for now.
}

void run_fit(vector< vector<double> > & data)
{
	// Setup the interface to multinest, run it.
	// Print out what parameters are being used here:
	printf("Starting fit with the following limits: \n");
	printf("Param       : Min        Max\n");
	printf("Omega (deg) : %1.4e %1.4e \n", Omega_min * RAD_TO_DEG, Omega_max * RAD_TO_DEG);
	printf("inc   (deg) : %1.4e %1.4e \n", inc_min * RAD_TO_DEG, inc_max * RAD_TO_DEG);
	printf("omega (deg) : %1.4e %1.4e \n", omega_min * RAD_TO_DEG, omega_max * RAD_TO_DEG);
	printf("e           : %1.4e %1.4e \n", e_min, e_max);
	printf("T (time)    : %1.4e %1.4e \n", T_min, T_max);
	printf("tau (time)  : %1.4e %1.4e \n", tau_min, tau_max);

	if(fit_motion)
	{
		printf("x0 (arb u)  : %1.4e %1.4e \n", x_0_min, x_0_max);
		printf("y0 (arb u)  : %1.4e %1.4e \n", y_0_min, y_0_max);
		printf("mu_x (u/t)  : %1.4e %1.4e \n", mu_x_min, mu_x_max);
		printf("mu_y (u/t)  : %1.4e %1.4e \n", mu_y_min, mu_y_max);
		printf("pi (arb u)  : %1.4e %1.4e \n", pi_min, pi_max);
	}

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

	// Now compute the (sometimes partial) priors:
	prior_Omega = 1.0 / scale_Omega;
	prior_inc = 1.0 / scale_inc;
	prior_omega = 1.0 / scale_omega;
	prior_alpha = 1.0 / scale_alpha;
	prior_e = 1.0 / scale_e;
	prior_tau = 1.0 / log(tau_max / tau_min);
	prior_T = 1.0 / log(T_max / T_min);
	prior_x_0 = 1.0 / scale_x_0;
	prior_y_0 = 1.0 / scale_y_0;
	prior_mu_x = 1.0 / scale_mu_x;
	prior_mu_y = 1.0 / scale_mu_y;
	prior_pi = 1.0 / scale_pi;


    n_data = data.size();
    printf("Found %i data points.\n", n_data);
    int opt_params = 0;

    if(fit_motion)
    	opt_params += 5;

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
	if(Omega_min == 0 && Omega_max == TWO_PI)
		pWrap[0] = 1;

	if(inc_min == 0 && inc_max == TWO_PI)
		pWrap[1] = 1;

	if(omega_min == 0 && omega_max == TWO_PI)
		pWrap[2] = 1;


	char root[100] = "chains/fitast-";		// root for output files
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
		{
			printf("NOTE: Fitting zero points, proper motions, and parallax.\n");
			fit_motion = true;
		}

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

		if(strcmp(argv[i], "-err") == 0)
			default_error = atof(argv[i + 1]);

    }
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

    if(argc < 2)
    	cout << "Missing filename on command line";


    // Read in the input filename:
    string input_rv = string(argv[1]);

    // Parse remaining, common parameters.
    ParseCommonParams(argc, argv, omega_min, omega_max, e_min, e_max, tau_min, tau_max, T_min, T_max, param_error);
    ParseProgOptions(argc, argv, param_error);

    if(param_error)
    	return 0;

    // Read in the file of RV data.
    // Each row should contain time, RV, [errors]
    const string comment_chars("\\#~$&£%");

    vector< vector<int> > split_info;

    read_data(input_rv, comment_chars, default_error, split_info, data);

    run_fit(data);

	return 0;
}
