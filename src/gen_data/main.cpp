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
#include "orbit.h"
#include "common.h"

using namespace std;

// Prints out help describing the options on the command line
void PrintHelp()
{
	string usage = "gendata - RV/Astrometry Test Dataset Generator\n"
	"Creates random sample orbital data to test the fitting routines\n\n"
	"Usage: \n"
	" gendata output_file_rv output_file_astrometry output_file_parameters [optional parameters]\n"
	" \n"
	"Optional Arguments: \n"
	" -h     Prints this message \n"
	"";

	cout << usage << "\n";
}

// Generates some orbital data.
void GenerateData(double Omega, double inc, double omega, double a, double alpha, double e, double tau, double T,
		vector<double> & times, vector< vector<double> > & positions, vector<double> & velocities)
{
	int n_data = times.size();
	double t, x, y, z, rv;

	for(int i = 0; i < n_data; i++)
	{
		// Compute the positions and velocities
		t = times[i];
		GetPositions(Omega, inc, omega, alpha, e, tau, T, t, x, y, z);
		GetRV(inc, omega, a, e, tau, T, t, rv);

		positions[i][0] = x;
		positions[i][1] = y;
		velocities[i] = rv;
	}
}

// Writes the data and parameter information to files.
void WriteData(double Omega, double inc, double omega, double a, double alpha, double e, double tau, double T,
		vector<double> times, vector< vector<double> > positions, vector<double> velocities,
		string output_rv, string output_astr, string output_params)
{
	// First write out the parameter file:
	ofstream params;
	params.precision(10);
	params.setf(ios::fixed,ios::floatfield);
	params.open(output_params.c_str());
	params << "Omega = " << Omega << endl;
	params << "inc   = " << inc << endl;
	params << "omega = " << omega << endl;
	params << "a     = " << a << endl;
	params << "alpha = " << alpha << endl;
	params << "e     = " << e << endl;
	params << "tau   = " << tau << endl;
	params << "T     = " << T << endl;
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

	int n_data = times.size();
	double t;
	for(int i = 0; i < n_data; i++)
	{
		t = times[i];
		// Write out the RV value
		rv << t << " " << velocities[i] << endl;
		astr << t << " " << RadToMas(positions[i][0]) << " " << RadToMas(positions[i][1]) << endl;
	}

	// Now close the files.
	rv.close();
	astr.close();


}

// The main routine.  Basically just used to parse out some parameters before handing
// things off to other functions.
int main(int argc, char *argv[])
{
    if(argc == 1)
    {
        PrintHelp();
        return 0;
    }

    if(argc < 4)
    	cout << "Missing filename on command line";

    // Read in the output filenames:
    string output_rv = string(argv[1]);
    string output_ast = string(argv[2]);
    string output_param = string(argv[3]);

    // Now generate some random orbital parameters.  Angular units in degrees for now.
    double Omega = 45;
    double inc = 20;
    double omega = 90;
    double a = 5.2E3;
    double alpha = 30;
    double e = 0.227;
    double tau = 580;
    double T = 9890;

    // How many data points?
    int n_data = 100;
    double t_0 = 0;
    double dt = 100;

	for (int i = 3; i < argc; i++)
	{
		// First see if the user is requesting help:
		if(strcmp(argv[i], "-h") == 0)
		{
			PrintHelp();
		}

// This is left here just in case we need a template for parsing any additional parameters
//		// We need to know some information about the target:
//		if ((strcmp(argv[i], "-t") == 0) && (i < argc - 1))
//		{
//			target.ImportFile(string(argv[i+1]), comment_chars);
//			target.ParseFileOptions(argv, i+2, argc);
//			n_params += 1;
//		}

	}

	// Now convert the angular orbital paramters over to radians
	Omega = DegToRad(Omega);
	inc = DegToRad(inc);
	omega = DegToRad(omega);
	alpha = MasToRad(alpha);

	// Create vectors into which the data will be stored.
	vector< double > times;
	vector< vector<double> > positions;
	vector< double > velocities;

	// Resize the vectors:
	times.resize(n_data);
	positions.resize(n_data);
	velocities.resize(n_data);
	for(int i = 0; i < n_data; i++)
	{
		// Calculate the observation time
		times[i] = t_0 + i * dt;

		// Resize the position vector, store only (x,y) pairs.
		positions[i].resize(2);
	}

	// Now generate the data and write it out to a file.
	GenerateData(Omega, inc, omega, a, alpha, e, tau, T, times, positions, velocities);
	WriteData(Omega, inc, omega, a, alpha, e, tau, T, times, positions, velocities, output_rv, output_ast, output_param);

	return 0;
}
