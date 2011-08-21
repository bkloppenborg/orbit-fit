/*
 * main.cpp
 *
 *  Created on: Aug 14, 2011
 *      Author: bkloppenborg
 *
 *  Generate random sample data to test the fitting algorithms.
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

#include "multinest.h"

using namespace std;

// Global Variables (scales and (sometimes partial) priors:
double V0_min;
double V0_max;
double omega_min;
double omega_max;
double asini_min;
double asini_max;
double e_min;
double e_max;
double tau_min;
double tau_max;
double T_min;
double T_max;
double s_min;
double s_max;

double scale_V0;
double scale_omega;
double scale_asini;
double scale_e;
double scale_tau;
double scale_T;
double scale_s;

double prior_V0;
double prior_omega;
double prior_asini;
double prior_e;
double prior_tau;
double prior_T;
double prior_s;

#define DAY_IN_SEC 24 * 60 * 60

vector< vector<double> > data;
int n_data;

// Prints out help describing the options on the command line
void print_help()
{
	string usage = "fitrv - A Bayesian Radial Velocity Fitting Program\n"
	"Fits Radial Velocity data using MultiNest\n\n"
	"Usage: \n"
	" fitrv input_file \n"
	" \n"
	"Optional Arguments: \n"
	" -h     Prints this message \n"
	"";

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
			// Covert rv over to km/day
			rv = atof(tokens[1].c_str()) * DAY_IN_SEC;
		}
		catch(...)
		{
			throw std::runtime_error("Could not parse line in RV data file.");
		}

		// Now do the RV, permit
		try
		{
			e_rv = atof(tokens[2].c_str()) * DAY_IN_SEC;
		}
		catch(...)
		{
			e_rv = 0;
		}

		if(e_rv == 0)
			e_rv = default_error * DAY_IN_SEC;

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
	double t, V, rvi, e_rvi, e_rvi2, tmp;

    // First extract the parameters, convert to real units:
	double V0 = Cube[0] * scale_V0 + V0_min;
    double omega = Cube[1] * scale_omega;
    double asini = Cube[2] * scale_asini;
    double e = Cube[3]; // not scaled
    double tau = Cube[4] * scale_tau;
    double T = Cube[5] * scale_T;
    double s = Cube[6] * scale_s;

    // Now set the cube parameters:
    Cube[0] = V0;
    Cube[1] = omega * RAD_TO_DEG;
    Cube[2] = asini;
    Cube[3] = e;
    Cube[4] = tau;
    Cube[5] = T;
    Cube[6] = s;

    // Compute a few things
    double prior = prior_V0
    			+ prior_omega
    			+ 1.0 / log(asini) * prior_asini
    			+ prior_e
    			+ 1.0 / log(tau) * prior_tau
    			+ 1.0 / log(T) * prior_T;

    double s2 = s*s; // intrinsic noise due to star or whatever...
    double llike = (double) n_data / 2.0 * log(TWO_PI);

    // Now compute the contribution of loglike from the data - model:
    for(int i = 0; i < n_data; i++)
    {
        t = data[i][0];
        rvi = data[i][1];
        e_rvi = data[i][2];

        GetRV(omega, asini, e, tau, T, t, V);
        e_rvi2 = e_rvi * e_rvi;

        tmp = e_rvi2 + s2;

        llike -= 0.5 * log(tmp) + (V0+V-rvi)*(V0+V-rvi) / (2*tmp);
    }

    //cout << omega << " " << asini << " " << e << " " << tau << " " << T << " " << llike << endl;

	// Assign the value and we're done.
	*lnew = llike + prior;
}

void run_fit(vector< vector<double> > & data)
{
	// Setup the interface to multinest, run it.
    // TODO: Read in these limits from elsewhere
	V0_min = -10;
	V0_max = 10;
    omega_min = 0;
    omega_max = TWO_PI;
    asini_min = 0;
    asini_max = 1E10;
    e_min = 0;
    e_max = 1;
    tau_min = 0;
    tau_max = 2.5E7;
    T_min = 0;
    T_max = 1E5;
    s_min = 0;
    s_max = 15;

    // Compute scales:
    scale_V0 = V0_max - V0_min;
    scale_omega = omega_max - omega_min;
    scale_asini = asini_max - asini_min;
    scale_e = e_max - e_min;
    scale_tau = tau_max - tau_min;
    scale_T = T_max - T_min;
    scale_s = s_max - s_min;

    // Compute much of values needed for the priors:
    // Here omega and e will use uniform priors
    // asini, tau, and T use Jeffrey's priors.
    prior_V0 = 1.0 / scale_V0;
    prior_omega = 1.0 / scale_omega;
    prior_asini = 1.0 / log(asini_max / asini_min);
    prior_e = 1.0 / scale_e;
    prior_tau = 1.0 / log(tau_max / tau_min);
    prior_T = 1.0 / log(T_max / T_min);
    prior_s = 1.0 / scale_s;

    n_data = data.size();
    printf("Found %i data points.\n", n_data);


	// set the MultiNest sampling parameters
	int mmodal = 1;					// do mode separation?
	int ceff = 0;					// run in constant efficiency mode?
	int nlive = 1000;				// number of live points
	double efr = 1.0;				// set the required efficiency
	double tol = 0.5;				// tol, defines the stopping criteria
	int ndims = 7;					// dimensionality (no. of free parameters)
	int nPar = 7;					// total no. of parameters including free & derived parameters
	int nClsPar = 7;				// no. of parameters to do mode separation on
	int updInt = 100;				// after how many iterations feedback is required & the output files should be updated
									// note: posterior files are updated & dumper routine is called after every updInt*10 iterations
	double Ztol = -1E90;			// all the modes with logZ < Ztol are ignored
	int maxModes = 100;				// expected max no. of modes (used only for memory allocation)
	int pWrap[ndims];				// which parameters to have periodic boundary conditions?
	for(int i = 0; i < ndims; i++)
	    pWrap[i] = 0;

	pWrap[1] = 1;   // omega

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

void dumper(int *nSamples, int *nlive, int *nPar, double **physLive,
    double **posterior, double *paramConstr, double *maxLogLike, double *logZ)
{
    // Do nothing.
}

// A wrapper to kick off multinest.
void run_multinest(int mmodal, int ceff, int nlive, double tol,
    double efr, int ndims, int nPar, int nClsPar,
    int maxModes, int updInt, double Ztol, char root[],
    int seed, int *pWrap, int fb, int resume,
    int outfile, int initMPI, double logZero,
    void (*LogLike)(double *, int *, int *, double *),
    void (*dumper)(int *, int *, int *, double **, double **, double *, double *, double *),
    int context)
{
    // Clear out the remaining characters in the string:
    int i;
	for (i = strlen(root); i < 100; i++)
	    root[i] = ' ';

    // Run the nested sampling algorithm
    NESTRUN(&mmodal, &ceff, &nlive, &tol,
        &efr, &ndims, &nPar, &nClsPar,
        &maxModes, &updInt, &Ztol, root,
        &seed, pWrap, &fb, &resume,
        &outfile, &initMPI, &logZero,
        LogLike,
        dumper,
        &context);
}

// The main routine.  Basically just used to parse out some parameters before handing
// things off to other functions.
int main(int argc, char *argv[])
{
    if(argc == 1)
    {
        print_help();
        return 0;
    }

    if(argc < 2)
    	cout << "Missing filename on command line";

    // Read in the input filename:
    string input_rv = string(argv[1]);

    // Read in the file of RV data.
    // Each row should contain time, RV, [errors]
    const string comment_chars("\\#~$&£%");

    int pos_size[] = {0, 12, 13, 6, 20, 4};
    vector< vector<int> > split_info;
    for(int i = 0; i < 3; i++)
    {
    	vector<int> tmp;
    	tmp.push_back(pos_size[2*i]);
    	tmp.push_back(pos_size[2*i+1]);
    	split_info.push_back(tmp);
    }

    read_data(input_rv, comment_chars, 5, split_info, data);

    run_fit(data);

	return 0;
}