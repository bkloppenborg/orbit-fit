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

// Importing/exporting unit to radian conversion factors
double UNIT_TO_RAD = 1;
double RAD_TO_UNIT = 1;

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
			t = atof(tokens.at(0).c_str()) * DAY_TO_SEC;
			x = atof(tokens.at(1).c_str()) * UNIT_TO_RAD;
			y = atof(tokens.at(3).c_str()) * UNIT_TO_RAD;
		}
		catch(...)
		{
			throw std::runtime_error("Could not parse line in astrometric data file.");
		}

		// Now for the uncertainties (these may not exist)
		try
		{
			e_x = atof(tokens.at(2).c_str()) * UNIT_TO_RAD;
		}
		catch(...)
		{
			e_x = 0;
		}

		try
		{
			e_y = atof(tokens.at(4).c_str()) * UNIT_TO_RAD;
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
			e_x = default_error * UNIT_TO_RAD;
		if(e_y == 0)
			e_y = default_error * UNIT_TO_RAD;

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
	// Local variables
	double t, xi, e_xi, yi, e_yi, P_alpha, P_delta;
	double x, y, z, err, dt;
	double e_xi2, e_yi2;

	// Used only if -motion is specified on the command line.
	double x_0, y_0, mu_x, mu_y, pi;

    // First extract the parameters, convert to real units:
    double Omega = 	Cube[0] * scale_Omega + Omega_min;
    double inc = 	Cube[1] * scale_inc + inc_min;
    double omega = 	Cube[2] * scale_omega + omega_min;
    double alpha =	Cube[3] * scale_alpha + alpha_min;
    double e = 		Cube[4] * scale_e + e_min;
    double tau = 	Cube[5] * scale_tau + tau_min;
    double T = 		Cube[6] * scale_T + T_min;

	if(fit_motion)
	{
		x_0 = Cube[7] * scale_x_0 + x_0_min;
		y_0 = Cube[8] * scale_y_0 + y_0_min;
		mu_x = Cube[9] * scale_mu_x + mu_x_min;
		mu_y = Cube[10] * scale_mu_y + mu_y_min;
		pi = Cube[11] * scale_pi + pi_min;
	}

    // Now set the cube parameters:
    Cube[0] = Omega * RAD_TO_DEG;
    Cube[1] = inc * RAD_TO_DEG;
    Cube[2] = omega * RAD_TO_DEG;
    Cube[3] = alpha * RAD_TO_MAS;
    Cube[4] = e;
    Cube[5] = tau * SEC_TO_DAY;
    Cube[6] = T * SEC_TO_DAY;

    if(fit_motion)
    {
		Cube[7] = x_0 * RAD_TO_UNIT;
		Cube[8] = y_0 * RAD_TO_UNIT;
		Cube[9] = mu_x * RADSEC_TO_MASYR;
		Cube[10] = mu_y * RADSEC_TO_MASYR;
		Cube[11] = pi * RAD_TO_MAS;
    }

    // Compute a few things
    double prior = prior_Omega
    			+ prior_inc
    			+ prior_omega
    			+ 1.0 / alpha * prior_alpha
    			+ prior_e
    			+ 1.0 / tau * prior_tau
    			+ 1.0 / T * prior_T;

    if(fit_motion)
    	prior += prior_x_0 + prior_y_0 + prior_mu_x + prior_mu_y + prior_pi;

    double llike = 0;

    // Now compute the contribution of loglike from the data - model:
    for(int i = 0; i < n_data; i++)
    {
        t = data[i][0];
        xi = data[i][1];
        e_xi = data[i][2];
        yi = data[i][3];
        e_yi = data[i][4];
        P_alpha = data[i][5];
        P_delta = data[i][6];

        e_xi2 = e_xi * e_xi;
        e_yi2 = e_yi * e_yi;

        // Get the positions, compute the residuals.
        GetPositions(Omega, inc, omega, alpha, e, tau, T, t, x, y, z);

        if(fit_motion)
        {
        	dt = t - tau;
        	x += x_0 + mu_x * dt + pi * P_alpha;
        	y += y_0 + mu_y * dt + pi * P_delta;

        	//cout << dt << " " << x << " " << y << endl;
        }

        llike -= 0.5 * log(TWO_PI * e_xi * e_yi) + (x - xi)*(x - xi) / (2 * e_xi2) + (y - yi)*(y - yi) / (2 * e_yi2);
    }

    //cout << Omega << " " << inc << " " << omega << " " << alpha << " " << e << " " << tau << " " << T << " " << llike << endl;
    //cout << "X: " << x - xi << " " << e_xi << " Y: " << y - yi << " " << e_yi << endl;

	// Assign the value and we're done.
	*lnew = llike + prior;
}

void run_fit(vector< vector<double> > & data)
{
	// Setup the interface to multinest, run it.

	// Print out what parameters are being used here:
	printf("Starting fit with the following limits: \n");
	printf("Param: Min Max\n");
	printf("Omega: %f %f \n", Omega_min * RAD_TO_DEG, Omega_max * RAD_TO_DEG);
	printf("inc:   %f %f \n", inc_min * RAD_TO_DEG, inc_max * RAD_TO_DEG);
	printf("omega: %f %f \n", omega_min * RAD_TO_DEG, omega_max * RAD_TO_DEG);
	printf("alpha: %f %f \n", alpha_min,alpha_max);
	printf("e:     %f %f \n", e_min, e_max);
	printf("tau:   %e %e \n", tau_min * SEC_TO_DAY, tau_max * SEC_TO_DAY);
	printf("T:     %e %e \n", T_min * SEC_TO_DAY, T_max * SEC_TO_DAY);

	// Now spit out the optional parameters, if used.
    if(fit_motion)
    {
    	printf("x_0:   %f %f \n", x_0_min * RAD_TO_UNIT, x_0_max * RAD_TO_UNIT);
    	printf("y_0:   %f %f \n", y_0_min * RAD_TO_UNIT, y_0_max * RAD_TO_UNIT);
    	printf("mu_x:  %f %f \n", mu_x_min * RADSEC_TO_MASYR, mu_x_max * RADSEC_TO_MASYR);
    	printf("mu_y:  %f %f \n", mu_y_min * RADSEC_TO_MASYR, mu_y_max * RADSEC_TO_MASYR);
    	printf("pi:    %f %f \n", pi_min * RAD_TO_MAS, pi_max * RAD_TO_MAS);
    }

    scale_x_0 = x_0_max - x_0_min;
    scale_y_0 = y_0_max - y_0_min;
    scale_mu_x = mu_x_max - mu_x_min;
    scale_mu_y = mu_y_max - mu_y_min;
    scale_pi = pi_max - pi_min;
	scale_Omega = Omega_max - Omega_min;
	scale_inc = inc_max - inc_min;
	scale_omega = omega_max - omega_min;
	scale_alpha = alpha_max - alpha_min;
	scale_e = e_max - e_min;
    scale_tau = tau_max - tau_min;
    scale_T = T_max - T_min;

    // Motion information all have uniform priors:
    prior_x_0 = 1.0 / scale_x_0;
    prior_y_0 = 1.0 / scale_y_0;
    prior_mu_x = 1.0 / scale_mu_x;
    prior_mu_y = 1.0 / scale_mu_y;
    prior_pi = 1.0 / scale_pi;

	prior_Omega = 1.0 / scale_Omega;
	prior_inc = 1.0 / scale_inc;
	prior_omega = 1.0 / scale_omega;
	prior_alpha = 1.0 / log(alpha_max / alpha_min);
	prior_e = 1.0 / scale_e;
	prior_tau = 1.0 / log(tau_max / tau_min);
	prior_T = 1.0 / log(T_max / T_min);

    n_data = data.size();
    printf("Found %i data points.\n", n_data);

    // Check the optional parameters:
    if(fit_motion)
    	opt_params += 5;

	// set the MultiNest sampling parameters
	int mmodal = 1;					// do mode separation?
	int ceff = 0;					// run in constant efficiency mode?
	int nlive = 1000;				// number of live points
	double efr = 1.0;				// set the required efficiency
	double tol = 0.5;				// tol, defines the stopping criteria
	int ndims = 7 + opt_params;					// dimensionality (no. of free parameters)
	int nPar = 7 + opt_params;					// total no. of parameters including free & derived parameters
	int nClsPar = 7 + opt_params;				// no. of parameters to do mode separation on
	int updInt = 100;				// after how many iterations feedback is required & the output files should be updated
									// note: posterior files are updated & dumper routine is called after every updInt*10 iterations
	double Ztol = -1E90;			// all the modes with logZ < Ztol are ignored
	int maxModes = 100;				// expected max no. of modes (used only for memory allocation)
	int pWrap[ndims];				// which parameters to have periodic boundary conditions?
	for(int i = 0; i < ndims; i++)
	    pWrap[i] = 0;

	pWrap[0] = 1;	// Omega
	pWrap[2] = 1;   // omega

	if(fit_motion)
	{
		pWrap[7] = 1;	// x_0
	}

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
	x_0_min = 0;
	x_0_max = TWO_PI;
	y_0_min = -90 * DEG_TO_RAD;
	y_0_max = 90 * DEG_TO_RAD;
	mu_x_min = 0;
	mu_x_max = 20 * MASYR_TO_RADSEC;
	mu_y_min = 0;
	mu_y_max = 20 * MASYR_TO_RADSEC;
	pi_min = 0;
	pi_max = 10 * MAS_TO_RAD;

	for (int i = 1; i < argc; i++)
	{
		// First see if the user is requesting help:
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
		{
			mu_x_min = atof(argv[i+1]) * MASYR_TO_RADSEC;
		}

		if(strcmp(argv[i], "-mu_x_max") == 0)
		{
			mu_x_max = atof(argv[i+1]) * MASYR_TO_RADSEC;
		}

		if(strcmp(argv[i], "-mu_y_min") == 0)
		{
			mu_y_min = atof(argv[i+1]) * MASYR_TO_RADSEC;
		}

		if(strcmp(argv[i], "-mu_y_max") == 0)
		{
			mu_y_max = atof(argv[i+1]) * MASYR_TO_RADSEC;
		}

		if(strcmp(argv[i], "-pi_min") == 0)
		{
			pi_min = atof(argv[i+1]) * MAS_TO_RAD;
		}

		if(strcmp(argv[i], "-pi_max") == 0)
		{
			pi_max = atof(argv[i+1]) * MAS_TO_RAD;
		}

		if(strcmp(argv[i], "-x0_min") == 0)
		{
			x_0_min = atof(argv[i+1]) * UNIT_TO_RAD;
		}

		if(strcmp(argv[i], "-x0_max") == 0)
		{
			x_0_max = atof(argv[i+1]) * UNIT_TO_RAD;
		}

		if(strcmp(argv[i], "-y0_min") == 0)
		{
			y_0_min = atof(argv[i+1]) * UNIT_TO_RAD;
		}

		if(strcmp(argv[i], "-y0_max") == 0)
		{
			y_0_max = atof(argv[i+1]) * UNIT_TO_RAD;
		}

		if(strcmp(argv[i], "-units") == 0)
		{
			// valid options: mas, deg, rad
			string tmp = string(argv[i+1]);

			if(tmp == "rad")
				UNIT_TO_RAD = 1;
			else if(tmp == "deg")
				UNIT_TO_RAD = DEG_TO_RAD;
			else if(tmp == "mas")
				UNIT_TO_RAD = MAS_TO_RAD;
			else
			{
				printf("Conversion of input units, %s, to radians is not implemented.  Please use rad, deg, or mas.", tmp.c_str());
				param_error = true;
			}

			RAD_TO_UNIT = 1.0 / UNIT_TO_RAD;
		}


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

	ParseProgOptions(argc, argv, param_error);

    // Read in the input filename:
    string input_rv = string(argv[1]);

    // Parse remaining, common parameters.
    ParseCommandLine(argc, argv, Omega_min, Omega_max, inc_min, inc_max, omega_min, omega_max, tmp, tmp,
    		alpha_min, alpha_max, e_min, e_max, tau_min, tau_max, T_min, T_max, param_error);

    if(param_error)
    	return 0;

    // Read in the file of RV data.
    // Each row should contain time, RV, [errors]
    const string comment_chars("\\#~$&£%");

//    int pos_size[] = {0, 12, 13, 6, 20, 4};
    vector< vector<int> > split_info;
//    for(int i = 0; i < 3; i++)
//    {
//    	vector<int> tmp;
//    	tmp.push_back(pos_size[2*i]);
//    	tmp.push_back(pos_size[2*i+1]);
//    	split_info.push_back(tmp);
//    }

    double default_error = MasToRad(5.0);

    read_data(input_rv, comment_chars, default_error, split_info, data);

    run_fit(data);

	return 0;
}
