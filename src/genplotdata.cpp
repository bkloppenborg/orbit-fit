/*
 * genplotdata.cpp
 *
 *  Created on: Sep 19, 2011
 *      Author: bkloppenborg
 */

#include "genplotdata.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "constants.h"
#include "orbit.h"
#include "ReadTextFile.h"
#include "common.h"
#include "read_data.h"

using namespace std;

vector<vector<double> > read_data_summary(string filename, string comment_chars)
{
	// Read in the file, return it as a vector of a vector of doubles.
    vector<string> lines = ReadFile(filename, comment_chars, "Could not open summary file.");
    vector<vector<double> > results;
	vector<double> temp;
	vector< string > tokens;

	for (unsigned int i = 0; i < lines.size(); i++)
	{
		// Empty the storage vectors
		temp.clear();
		tokens.clear();

		tokens = Tokenize(lines[i]);

		// Convert the numbers over to doubles and push them onto the results vector
		for(unsigned int j = 0; j < tokens.size(); j++)
			temp.push_back(atof(tokens[j].c_str()));

		results.push_back(temp);
	}

	return results;
}

void extract_params_ast(vector<double> params, bool & fit_pm,
		double & omega, double & e, double & tau, double & T,
		double & Omega, double & inc, double & alpha,
		double & x0, double & y0, double & mu_x, double & mu_y, double & pi)
{
	// First inspect the length of the parameter vector to see if proper motions were fit:
	// The *-summary.txt file has nPar*4+2 columns.  Find nPar
	int n_par = (params.size() - 2) / 4;

	if(n_par == 12)
		fit_pm = true;

	omega = params[0] * DEG_TO_RAD;
	e = params[1];
	tau = params[2];
	T = params[3];
	Omega = params[4] * DEG_TO_RAD;
	inc = params[5] * DEG_TO_RAD;
	alpha = params[6];

	if(fit_pm)
	{
		x0 = params[7];
		y0 = params[8];
		mu_x = params[9];
		mu_y = params[10];
		pi = params[11];
	}
}

void extract_params_rv(vector<double> params,
		double & omega, double & e, double & tau, double & T,
		double & K, double & gamma)
{
	// We always fit six parameters.  Sometimes noise is included, but it isn't an important parameter here.
	omega = params[0] * DEG_TO_RAD;
	e = params[1];
	tau = params[2];
	T = params[3];
	K = params[4];
	gamma = params[5];
}


void gendata_ast(string output_basename, vector<vector<double> > params, double t_min, double t_max, double t_step)
{
	// With the astrometry data we either fit 7, 8, 12, or 13 parameters in fitrv, but the noise parameter
	// doesn't change anything in this function.
	// Inspect the first element to see if we have 12 parameters => proper motion fitting happens.
	double omega, e, tau, T, Omega, inc, alpha;
	double x0, y0, mu_x, mu_y, pi;
	double n, M, E, cos_E, sin_E, beta;
	double l1, m1, l2, m2, n1, n2;
	double t, x, xi, y, yi;
	double dt;
	bool fit_pm = false;

	// Iterate over the possible modes stored in params and write out a unique residual data file for each mode.
	for(int mode = 0; mode < params.size(); mode++)
	{
		// Set the output filename, extract the parameters
		stringstream output_name;
		output_name << output_basename << "_pos_" << mode + 1;

		ofstream output;
		output.precision(10);
		output.setf(ios::fixed,ios::floatfield);
		output.open(output_name.str().c_str());

		extract_params_ast(params[mode], fit_pm, omega, e, tau, T, Omega, inc, alpha, x0, y0, mu_x, mu_y, pi);

	    n = ComputeN(T);
		Compute_Coefficients(Omega, inc, omega, l1, m1, n1, l2, m2, n2);

	    for(double t = t_min; t < t_max; t += t_step)
	    {
			M = ComputeM(tau, n, t);
			E = ComputeE(M, e);

			cos_E = cos(E);
			sin_E = sin(E);
			beta = sqrt(1 - e*e);

			Compute_xy(alpha, beta, e, l1, l2, m1, m2, cos_E, sin_E, x, y);

			if(fit_pm)
			{
	    		dt = t - tau;
	    		x += x0 + mu_x * dt; // + pi * P_a;
	    		y += y0 + mu_y * dt; //  + pi * P_d;
			}

			output << t << " " << x << " " << y << endl;
	    }

	    output.close();

	}
}

void gendata_rv(string output_basename, vector<vector<double> > params, double t_min, double t_max, double t_step)
{
	// With the astrometry data we either fit 7, 8, 12, or 13 parameters in fitrv, but the noise parameter
	// doesn't change anything in this function.
	// Inspect the first element to see if we have 12 parameters => proper motion fitting happens.
	double omega, e, tau, T, Omega, K, gamma;
	double n, M, E, cos_E, sin_E, beta;
	double sin_omega, cos_omega;
	double l1, m1, l2, m2, n1, n2;
	double t, rv, rvi;

	// Iterate over the possible modes stored in params and write out a unique residual data file for each mode.
	for(int mode = 0; mode < params.size(); mode++)
	{
		// Set the output filename, extract the parameters
		stringstream output_name;
		output_name << output_basename << "_rv_" << mode + 1;

		ofstream output;
		output.precision(10);
		output.setf(ios::fixed,ios::floatfield);
		output.open(output_name.str().c_str());

		extract_params_rv(params[mode], omega, e, tau, T, K, gamma);

	    n = ComputeN(T);
	    sin_omega = sin(omega);
	    cos_omega = cos(omega);

	    for(double t = t_min; t < t_max; t += t_step)
	    {
			M = ComputeM(tau, n, t);
			E = ComputeE(M, e);

			cos_E = cos(E);
			sin_E = sin(E);
			beta = sqrt(1 - e*e);

	    	rv = K * beta / (1 - e * cos_E) * (beta * cos_E * cos_omega - sin_E * sin_omega);

			output << t << " " << rv << endl;
	    }

	    output.close();

	}
}

int main(int argc, char *argv[])
{
	bool param_error = false;
	bool astrometry = false;
	bool rv = false;

	if(argc == 1)
		print_help();

	if(argc < 5)
		param_error = true;

	if(param_error)
		return 0;

    // Read in the input filename:
    string output_type = string(argv[1]);
    string summary_file = string(argv[2]);
    string data_file = string(argv[3]);
    string output_basename = string(argv[4]);

    // Check that one of the programs was defined:
    if(output_type == "fitast")
    	astrometry = true;
    else if(output_type == "fitrv")
    	rv = true;
    else
    {
    	cout << "Data type must be either 'fitast' or 'fitrv'" << endl;
    	return 0;
    }

    // Read in the parameters
    const string comment_chars("\\#~$&£%");
    vector< vector<double> > params = read_data_summary(summary_file, comment_chars);

    // Now read in the data
    vector< vector<int> > split_info;	// intentionally blank, could be used for column-formatted data with missing entries.
    vector< vector<double> > data;

    // TODO: Fix the read_no_error, default_error, read_r_theta parameters here:

    if(astrometry)
    	read_data_ast(data_file, comment_chars, split_info, data, false, 0);
    else if(rv)
    	read_data_rv(data_file, comment_chars, split_info, data, false, 0);

    int n_data = 1000;	// how many data points do we generate over the interval t_max:t_min?
    double t_min = data[0][0];
    double t_max = data[data.size() - 1][0];
    double t_step = (t_max - t_min) / n_data;

    // Now generate and write out the residuals
    if(astrometry)
    	residuals_ast(output_basename, params, data);
    if(rv)
    	residuals_rv(output_basename, params, data);

    // Lastly generate some overlapping data for plotting:
    if(astrometry)
    	gendata_ast(output_basename, params, t_min, t_max, t_step);
    if(rv)
    	gendata_rv(output_basename, params, t_min, t_max, t_step);

    return 0;
}


void print_help()
{
	string usage = "genplotdata - A program for making plot data\n"
	"Parses the output of multinest and generates data and residuals for plotting\n\n"
	"Usage: \n"
	" genplotdata fitast|fitrv summary_file data_file output_basename [options] ... \n"
	"";

	cout << usage << "\n";
}

void residuals_ast(string output_basename, vector<vector<double> > params, vector<vector<double> > data)
{
	// With the astrometry data we either fit 7, 8, 12, or 13 parameters in fitrv, but the noise parameter
	// doesn't change anything in this function.
	// Inspect the first element to see if we have 12 parameters => proper motion fitting happens.
	double omega, e, tau, T, Omega, inc, alpha;
	double x0, y0, mu_x, mu_y, pi;
	double n, M, E, cos_E, sin_E, beta;
	double l1, m1, l2, m2, n1, n2;
	double t, x, xi, y, yi;
	double dt;
	bool fit_pm = false;

	// Iterate over the possible modes stored in params and write out a unique residual data file for each mode.
	for(int mode = 0; mode < params.size(); mode++)
	{
		// Set the output filename, extract the parameters
		stringstream output_name;
		output_name << output_basename << "_res_" << mode + 1;

		ofstream output;
		output.precision(10);
		output.setf(ios::fixed,ios::floatfield);
		output.open(output_name.str().c_str());

		extract_params_ast(params[mode], fit_pm, omega, e, tau, T, Omega, inc, alpha, x0, y0, mu_x, mu_y, pi);

	    n = ComputeN(T);
		Compute_Coefficients(Omega, inc, omega, l1, m1, n1, l2, m2, n2);

	    for(int i = 0; i < data.size(); i++)
	    {
	    	t = data[i][0];
	    	xi = data[i][1];
	    	yi = data[i][3];

			M = ComputeM(tau, n, t);
			E = ComputeE(M, e);

			cos_E = cos(E);
			sin_E = sin(E);
			beta = sqrt(1 - e*e);

			Compute_xy(alpha, beta, e, l1, l2, m1, m2, cos_E, sin_E, x, y);

			if(fit_pm)
			{
	    		dt = t - tau;
	    		x += x0 + mu_x * dt; // + pi * P_a;
	    		y += y0 + mu_y * dt; //  + pi * P_d;
			}

			output << t << " " << x - xi << " " << y - yi << endl;
	    }

	    output.close();

	}
}

void residuals_rv(string output_basename, vector<vector<double> > params, vector<vector<double> > data)
{
	// With the astrometry data we either fit 7, 8, 12, or 13 parameters in fitrv, but the noise parameter
	// doesn't change anything in this function.
	// Inspect the first element to see if we have 12 parameters => proper motion fitting happens.
	double omega, e, tau, T, K, gamma;
	double n, M, E, cos_E, sin_E, beta;
	double sin_omega, cos_omega;
	double l1, m1, l2, m2, n1, n2;
	double t, rv, rvi;

	// Iterate over the possible modes stored in params and write out a unique residual data file for each mode.
	for(int mode = 0; mode < params.size(); mode++)
	{
		// Set the output filename, extract the parameters
		stringstream output_name;
		output_name << output_basename << "_res_" << mode + 1;

		ofstream output;
		output.precision(10);
		output.setf(ios::fixed,ios::floatfield);
		output.open(output_name.str().c_str());

		extract_params_rv(params[mode], omega, e, tau, T, K, gamma);

	    n = ComputeN(T);
	    sin_omega = sin(omega);
	    cos_omega = cos(omega);

	    for(int i = 0; i < data.size(); i++)
	    {
	    	t = data[i][0];
	    	rvi = data[i][1];

			M = ComputeM(tau, n, t);
			E = ComputeE(M, e);

			cos_E = cos(E);
			sin_E = sin(E);
			beta = sqrt(1 - e*e);

			rv = K * beta / (1 - e * cos_E) * (beta * cos_E * cos_omega - sin_E * sin_omega);

			output << t << " " << rv - rvi << endl;
	    }

	    output.close();

	}
}

