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
#include <cstdlib>
#include <fstream>
#include <math.h>
#include "orbit.h"
#include "common.h"
#include "random.h"
#include "constants.h"

using namespace std;

// Prints out help describing the options on the command line
void print_help()
{
	string usage = "gendata - RV/Astrometry Test Dataset Generator\n"
	"Creates random sample orbital data to test the fitting routines\n\n"
	"Usage: \n"
	" gendata base_file_name [optional parameters]\n"
	" \n"
	"Optional Arguments: \n"
	" -h     Prints this message \n"
	"";

	cout << usage << "\n";
}

// Parse command-line options that are specific to this program
void ParseProgOptions(int argc, char *argv[],
		double & Omega, double & inc, double & omega,
		double & alpha, double&  K,
		double & e, double & T, double & tau,
		double & x0, double & y0, double & mu_x, double & mu_y, double & pi,
		bool & param_error)
{
	// Init the random number generator.
	// Comment the RanInit line to make the code generate the same orbital parameters.
	static Rand_t random_seed;
	//RanInit(random_seed);

	Omega = Randouble(random_seed) * 360;
	inc = Randouble(random_seed) * 360 - 180;
	omega = Randouble(random_seed) * 360;
	alpha = Randouble(random_seed) * 100;
	K = Randouble(random_seed) * 30;
	e = Randouble(random_seed);
	T = Randouble(random_seed) * 20;
	tau = Randouble(random_seed) * 2000;

	x0 = Randouble(random_seed) * 360;
	y0 = Randouble(random_seed) * 180 - 90;
	mu_x = Randouble(random_seed) * 2E-5 - 1E5;
	mu_y = Randouble(random_seed) * 2E-5 - 1E5;
	pi = Randouble(random_seed) * 1E-5;

	for (int i = 1; i < argc; i++)
	{
		// Help
		if(strcmp(argv[i], "-h") == 0)
		{
			print_help();
			param_error = true;
		}

		// First see if the user is requesting help:
		if(strcmp(argv[i], "-Omega") == 0)
			Omega = atof(argv[i+1]);

		if(strcmp(argv[i], "-inc") == 0)
			inc = atof(argv[i+1]);

		if(strcmp(argv[i], "-omega") == 0)
			omega = atof(argv[i+1]);

		if(strcmp(argv[i], "-alpha") == 0)
			alpha = atof(argv[i+1]);

		if(strcmp(argv[i], "-K") == 0)
			K = atof(argv[i+1]);

		if(strcmp(argv[i], "-e") == 0)
			e = atof(argv[i+1]);

		if(strcmp(argv[i], "-T") == 0)
			T = atof(argv[i+1]);

		if(strcmp(argv[i], "-tau") == 0)
			tau = atof(argv[i+1]);

		if(strcmp(argv[i], "-x0") == 0)
			x0 = atof(argv[i+1]);

		if(strcmp(argv[i], "-y0") == 0)
			y0 = atof(argv[i+1]);

		if(strcmp(argv[i], "-mu_x") == 0)
			mu_x = atof(argv[i+1]);

		if(strcmp(argv[i], "-mu_y") == 0)
			mu_y = atof(argv[i+1]);

		if(strcmp(argv[i], "-pi") == 0)
			pi = atof(argv[i+1]);

    }

	// Convert angles over to radians
	Omega = Omega * DEG_TO_RAD;
	inc = inc * DEG_TO_RAD;
	omega = omega * DEG_TO_RAD;
}

// Generates some orbital data.
void GenerateData(double Omega, double inc, double omega, double K, double alpha, double e, double tau, double T,
		double x0, double y0, double mu_x, double mu_y, double pi,
		vector<double> & times, vector< vector<double> > & positions, vector< vector<double> > & positions_with_motion, vector< vector<double> > & velocities)
{
	int n_data = times.size();
	double t, x, s_x, y, s_y, z, rv, s_rv, sig, dt;

	// Pull up the random number generator.
	static Rand_t random_seed;

	for(int i = 0; i < n_data; i++)
	{
		// Compute the positions and velocities
		t = times[i];
		dt = t - times[0];
		GetPositions(Omega, inc, omega, alpha, e, tau, T, t, x, y, z);

		// Get the RV.  Units for time don't matter here.
		GetRV_K(K, omega, e, tau, T, t, rv);

		// TODO: This really isn't a good way of creating
		// data with noise as we're assuming 6 mas error on positioning, and 0.1 km/s on RV
		sig = Rangauss(random_seed);
		s_x = 0.1 * y * sig;
		s_y = 0.1 * x * sig;
		s_rv = 0.1 * rv * sig;

		positions[i][0] = x; // + s_x * Rangauss(random_seed);
		positions[i][1] = fabs(s_x);
		positions[i][2] = y; // + s_y * Rangauss(random_seed);
		positions[i][3] = fabs(s_y);

		// Positions with motion.  Note, parallax isn't considered yet.
		positions_with_motion[i][0] = x + x0 + mu_x * dt;
		positions_with_motion[i][1] = fabs(s_x);
		positions_with_motion[i][2] = y + y0 + mu_y * dt;
		positions_with_motion[i][3] = fabs(s_y);

		velocities[i][0] = rv + s_rv * Rangauss(random_seed);
		velocities[i][1] = fabs(s_rv);
	}
}

// Writes the data and parameter information to files.
void WriteData(double Omega, double inc, double omega, double K, double alpha, double e, double tau, double T,
		double x0, double y0, double mu_x, double mu_y, double pi,
		vector<double> times,
		vector< vector<double> > positions, vector< vector<double> > position_w_motion, vector< vector<double> > velocities,
		string output_rv, string output_astr, string output_astr_pm, string output_params)
{
	// First write out the parameter file:
	ofstream params;
	params.precision(10);
	params.setf(ios::fixed,ios::floatfield);
	params.open(output_params.c_str());
	params << "Omega = " << RadToDeg(Omega) << endl;
	params << "inc   = " << RadToDeg(inc) << endl;
	params << "omega = " << RadToDeg(omega) << endl;
	params << "K     = " << K << endl;
	params << "alpha = " << alpha << endl;
	params << "e     = " << e << endl;
	params << "tau   = " << tau << endl;
	params << "T     = " << T << endl;
	params << "x0    = " << x0 << endl;
	params << "y0    = " << y0 << endl;
	params << "mu_x  = " << mu_x << endl;
	params << "mu_y  = " << mu_y << endl;
	//params << "pi    = " << pi << endl;
	params.close();

	// Now write out the RV and astrometric data
	ofstream rv;
	rv.precision(10);
	rv.setf(ios::fixed,ios::floatfield);
	rv.open(output_rv.c_str());

	ofstream astr;
	astr.precision(10);
	astr.setf(ios::fixed,ios::floatfield);
	astr.open(output_astr.c_str());

	ofstream astr_pm;
	astr_pm.precision(10);
	astr_pm.setf(ios::fixed,ios::floatfield);
	astr_pm.open(output_astr_pm.c_str());

	int n_data = times.size();
	double t;
	for(int i = 0; i < n_data; i++)
	{
		t = times[i];
		// Write out the RV value
		rv << t << " " << velocities[i][0] << " " << velocities[i][1] << endl;
		astr << t << " " << positions[i][0] << " " << positions[i][1] << " " << positions[i][2] << " " << positions[i][3] << endl;
		astr_pm << t << " " << position_w_motion[i][0] << " " << position_w_motion[i][1] << " " << position_w_motion[i][2] << " " << position_w_motion[i][3] << endl;
	}

	// Now close the files.
	rv.close();
	astr.close();
	astr_pm.close();


}

// The main routine.  Basically just used to parse out some parameters before handing
// things off to other functions.
int main(int argc, char *argv[])
{
	bool param_error = false;

    if(argc < 2)
    {
    	cout << "Missing filename on command line";
    	param_error = true;
    }

    // Read in the output filenames:
    string base_name = string(argv[1]);
    string output_rv = base_name + "_rv";
    string output_ast = base_name + "_ast";
    string output_ast_pm = base_name + "_ast_pm";
    string output_param = base_name + "_params";

    // Now generate some random orbital parameters.  Angular units in degrees for now.
    double Omega = 0;
    double inc = 0;
    double omega = 0;
    double K = 0;
    double alpha = 0;
    double e = 0;
    double tau = 0;
    double T = 0;
    double x0 = 0;
    double y0 = 0;
    double mu_x = 0;
    double mu_y = 0;
    double pi = 0;
    double t_start;
    double t_end;

    ParseProgOptions(argc, argv, Omega, inc, omega, alpha, K, e, T, tau, x0, y0, mu_x, mu_y, pi, param_error);

    if(param_error)
    	return 0;

    // How many data points?
    int n_data = 100;
    double t_min = tau - 2*T;
    double t_max = tau + 2*T;
    double dt = (t_max - t_min) / n_data;

	// Create vectors into which the data will be stored.
	vector< double > times;
	vector< vector<double> > positions;
	vector< vector<double> > positions_w_motion;;
	vector< vector<double> > velocities;

	// Resize the vectors:
	times.resize(n_data);
	positions.resize(n_data);
	positions_w_motion.resize(n_data);
	velocities.resize(n_data);
	for(int i = 0; i < n_data; i++)
	{
		// Calculate the observation time
		times[i] = t_min + i * dt;

		// Resize the position vector, store only (x,y) pairs.
		positions[i].resize(4);
		positions_w_motion[i].resize(4);
		velocities[i].resize(2);
	}

	// Now generate the data and write it out to a file.
	GenerateData(Omega, inc, omega, K, alpha, e, tau, T, x0, y0, mu_x, mu_y, pi, times, positions, positions_w_motion, velocities);
	WriteData(Omega, inc, omega, K, alpha, e, tau, T, x0, y0, mu_x, mu_y, pi, times, positions, positions_w_motion, velocities, output_rv, output_ast, output_ast_pm, output_param);

	return 0;
}
