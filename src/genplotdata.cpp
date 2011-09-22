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

void extract_params_ast(vector<double> & params, bool & fit_pm,
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

void extract_params_both(vector<double> & params, bool & fit_pm,
		double & omega, double & e, double & tau, double & T,
		double & K, double & gamma,
		double & Omega, double & inc, double & alpha,
		double & x0, double & y0, double & mu_x, double & mu_y, double & pi)
{
	// First inspect the length of the parameter vector to see if proper motions were fit:
	// The *-summary.txt file has nPar*4+2 columns.  Find nPar
	int n_par = (params.size() - 2) / 4;
	int i = 8; // counter for offsets for additional parameters

	bool fit_rv_noise = false;

	if(n_par == 10)
		fit_rv_noise = true;

	if(n_par == 14)
		fit_pm = true;

	if(n_par == 15)
	{
		fit_rv_noise = true;
		fit_pm = true;
	}

	omega = params[0] * DEG_TO_RAD;
	e = params[1];
	tau = params[2];
	T = params[3];
	K = params[4];
	gamma = params[5];
	Omega = params[6] * DEG_TO_RAD;
	inc = params[7] * DEG_TO_RAD;
	alpha = params[8];

	if(fit_rv_noise)
		i++;

	if(fit_pm)
	{
		x0	 = params[i + 1];
		y0 	 = params[i + 2];
		mu_x = params[i + 3];
		mu_y = params[i + 4];
		pi   = params[i + 5];
	}
}

void extract_params_rv(vector<double> & params,
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


void gendata_ast(string output_basename, vector<vector<double> > & params, double t_min, double t_max, double t_step)
{
	// With the astrometry data we either fit 7, 8, 12, or 13 parameters in fitrv, but the noise parameter
	// doesn't change anything in this function.
	// Inspect the first element to see if we have 12 parameters => proper motion fitting happens.
	double omega, e, tau, T, Omega, inc, alpha;
	double x0, y0, mu_x, mu_y, pi;
	bool fit_pm = false;

	// Iterate over the possible modes stored in params and write out a unique residual data file for each mode.
	for(int mode = 0; mode < params.size(); mode++)
	{
		// Set the output filename, extract the parameters
		stringstream output_name;
		output_name << output_basename << "_" << mode + 1 << "_pos_ast" ;


		extract_params_ast(params[mode], fit_pm, omega, e, tau, T, Omega, inc, alpha, x0, y0, mu_x, mu_y, pi);

		gendata_ast(output_name.str(), fit_pm, omega, e, tau, T, Omega, inc, alpha, x0, y0, mu_x, mu_y, pi, t_min, t_max, t_step);
	}
}

void gendata_ast(string filename, bool fit_pm,
		double omega, double e, double tau, double T,
		double Omega, double inc, double alpha,
		double x0, double y0, double mu_x, double mu_y, double pi,
		double t_min, double t_max, double t_step)
{
	double n, M, E, cos_E, sin_E, beta;
	double l1, m1, l2, m2, n1, n2;
	double t, x, xi, y, yi;
	double dt;

	stringstream output_name;
	output_name << filename << "_orbit";

	ofstream orbit_out;
	orbit_out.precision(10);
	orbit_out.setf(ios::fixed,ios::floatfield);
	orbit_out.open(output_name.str().c_str());

	ofstream orbit_pm_out;

	if(fit_pm)
	{
		output_name.clear();//clear any bits set
		output_name.str(string());
		output_name << filename << "_pm";

		orbit_pm_out.precision(10);
		orbit_pm_out.setf(ios::fixed,ios::floatfield);
		orbit_pm_out.open(output_name.str().c_str());
	}

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

		orbit_out << t << " " << x << " " << y << endl;

		if(fit_pm)
		{
			dt = t - tau;
			x += x0 + mu_x * dt; // + pi * P_a;
			y += y0 + mu_y * dt; //  + pi * P_d;

			orbit_pm_out << t << " " << x << " " << y << endl;
		}


	}

	orbit_out.close();

	if(fit_pm)
		orbit_pm_out.close();
}

void gendata_both(string output_basename, vector<vector<double> > & params, double t_min, double t_max, double t_step)
{
	// With the astrometry data we either fit 7, 8, 12, or 13 parameters in fitrv, but the noise parameter
	// doesn't change anything in this function.
	// Inspect the first element to see if we have 12 parameters => proper motion fitting happens.
	double omega, e, tau, T, Omega, inc, alpha;
	double K, gamma;
	double x0, y0, mu_x, mu_y, pi;
	bool fit_pm = false;

	// Iterate over the possible modes stored in params and write out a unique residual data file for each mode.
	for(int mode = 0; mode < params.size(); mode++)
	{
		// Set the output filename, extract the parameters
		// Set the output filename, extract the parameters
		stringstream ast_output;
		ast_output << output_basename << "_"  << mode + 1 << "_pos_ast";

		stringstream rv_output;
		rv_output << output_basename << "_"  << mode + 1 << "_pos_rv";

		// Extract the parameters, then generate the data
		extract_params_both(params[mode], fit_pm, omega, e, tau, T, K, gamma, Omega, inc, alpha, x0, y0, mu_x, mu_y, pi);

		gendata_ast(ast_output.str(), fit_pm, omega, e, tau, T, Omega, inc, alpha, x0, y0, mu_x, mu_y, pi, t_min, t_max, t_step);
		gendata_rv(rv_output.str(), omega, e, tau, T, K, gamma, t_min, t_max, t_step);

	}
}

void gendata_rv(string output_basename, vector<vector<double> > & params, double t_min, double t_max, double t_step)
{
	// With the astrometry data we either fit 7, 8, 12, or 13 parameters in fitrv, but the noise parameter
	// doesn't change anything in this function.
	// Inspect the first element to see if we have 12 parameters => proper motion fitting happens.
	double omega, e, tau, T, Omega, K, gamma;

	// Iterate over the possible modes stored in params and write out a unique residual data file for each mode.
	for(int mode = 0; mode < params.size(); mode++)
	{
		// Set the output filename, extract the parameters
		stringstream output_name;
		output_name << output_basename << "_"  << mode + 1 << "_pos_rv" ;

		extract_params_rv(params[mode], omega, e, tau, T, K, gamma);

		gendata_rv(output_name.str(), omega, e, tau, T, K, gamma, t_min, t_max, t_step);
	}
}


void gendata_rv(string output_name,
		double omega, double e, double tau, double T,
		double K, double gamma,
		double t_min, double t_max, double t_step)
{

	double n, M, E, cos_E, sin_E, beta;
	double sin_omega, cos_omega;
	double l1, m1, l2, m2, n1, n2;
	double t, rv, rvi;

	ofstream output;
	output.precision(10);
	output.setf(ios::fixed,ios::floatfield);
	output.open(output_name.c_str());

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

int main(int argc, char *argv[])
{
	bool param_error = false;
	bool astrometry = false;
	bool both = false;
	bool rv = false;

	if(argc == 1)
		print_help();

	if(argc < 5)
		param_error = true;

	if(param_error)
		return 0;

    // Read in the input filename:
    string output_type = string(argv[1]);

    // Check that one of the programs was defined:
    if(output_type == "fitast")
    	astrometry = true;
    else if(output_type == "fitrv")
    	rv = true;
    else if(output_type == "fitboth")
    {
    	both = true;
    	rv = true;
    	astrometry = true;
    }
    else
    {
    	cout << "Data type must be either 'fitast', 'fitrv', or 'fitboth'" << endl;
    	return 0;
    }

    string summary_file = string(argv[2]);

    string rv_datafile;
    string ast_datafile;
    string output_basename;

    if(both)
    {
    	rv_datafile = string(argv[3]);
    	ast_datafile = string(argv[4]);
    	output_basename = string(argv[5]);
    }
    else if(rv)
    {
    	rv_datafile = string(argv[3]);
    	output_basename = string(argv[4]);
    }
    else if(astrometry)
    {
    	ast_datafile = string(argv[3]);
    	output_basename = string(argv[4]);
    }



    // Read in the parameters
    const string comment_chars("\\#~$&£%");
    vector< vector<double> > params = read_data_summary(summary_file, comment_chars);

    // Now read in the data
    vector< vector<int> > split_info;	// intentionally blank, could be used for column-formatted data with missing entries.
    vector< vector<double> > rv_data;
    vector< vector<double> > ast_data;

    // TODO: Fix the read_no_error, default_error, read_r_theta parameters here:

    if(astrometry)
    	read_data_ast(ast_datafile, comment_chars, split_info, ast_data, false, 0);
    if(rv)
    	read_data_rv(rv_datafile, comment_chars, split_info, rv_data, false, 0);

    // Now determine the min/max dates in the data sets:
    int n_data = 1000;
    double t_min, t_max, t_step;

    if(both)
    {
    	t_min = min(rv_data[0][0], ast_data[0][0]);
    	t_max = max(rv_data[rv_data.size() - 1][0], ast_data[rv_data.size() - 1][0]);
    }
    else if(astrometry)
    {
    	t_min = ast_data[0][0];
    	t_max = ast_data[ast_data.size() - 1][0];
    }
    else if(rv)
    {
    	t_min = rv_data[0][0];
    	t_max = rv_data[rv_data.size() - 1][0];
    }

    t_step = (t_max - t_min) / n_data;

    // Now generate and write out the residuals
    if(both)
    	residuals_both(output_basename, params, ast_data, rv_data);
    else if(astrometry)
    	residuals_ast(output_basename, params, ast_data);
    else if(rv)
    	residuals_rv(output_basename, params, rv_data);

    // Lastly generate some overlapping data for plotting:
    if(both)
    	gendata_both(output_basename, params, t_min, t_max, t_step);
    else if(astrometry)
    	gendata_ast(output_basename, params, t_min, t_max, t_step);
    else if(rv)
    	gendata_rv(output_basename, params, t_min, t_max, t_step);

    return 0;
}


void print_help()
{
	string usage = "genplotdata - A program for making plot data\n"
	"Parses the output of multinest and generates data and residuals for plotting\n\n"
	"Usage: \n"
	" genplotdata fitast|fitrv summary_file data_file output_basename [options] ... \n"
	"or\n"
	" genplotdata fitboth summary_file rv_data_file ast_data_file output_basename [options] ... \n"
	"\n\n";

	cout << usage << "\n";
}

void residuals_ast(string output_basename, vector<vector<double> > & params, vector<vector<double> > & data)
{
	// With the astrometry data we either fit 7, 8, 12, or 13 parameters in fitrv, but the noise parameter
	// doesn't change anything in this function.
	// Inspect the first element to see if we have 12 parameters => proper motion fitting happens.
	double omega, e, tau, T, Omega, inc, alpha;
	double x0, y0, mu_x, mu_y, pi;
	bool fit_pm = false;

	// Iterate over the possible modes stored in params and write out a unique residual data file for each mode.
	for(int mode = 0; mode < params.size(); mode++)
	{
		// Set the output filename, extract the parameters
		stringstream output_name;
		output_name << output_basename << "_"  << mode + 1 << "_res_ast";

		extract_params_ast(params[mode], fit_pm, omega, e, tau, T, Omega, inc, alpha, x0, y0, mu_x, mu_y, pi);

		residuals_ast(output_name.str(), fit_pm, omega, e, tau, T, Omega, inc, alpha, x0, y0, mu_x, mu_y, pi, data);
	}
}

void residuals_ast(string output_name, bool fit_pm,
		double omega, double e, double tau, double T,
		double Omega, double inc, double alpha,
		double x0, double y0, double mu_x, double mu_y, double pi,
		vector<vector<double> > & data)
{
	double n, M, E, cos_E, sin_E, beta;
	double l1, m1, l2, m2, n1, n2;
	double t, x, xi, y, yi;
	double dt;

	ofstream output;
	output.precision(10);
	output.setf(ios::fixed,ios::floatfield);
	output.open(output_name.c_str());

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

void residuals_both(string output_basename, vector<vector<double> > & params, vector<vector<double> > & ast_data, vector<vector<double> > & rv_data)
{
	// With the astrometry data we either fit 7, 8, 12, or 13 parameters in fitrv, but the noise parameter
	// doesn't change anything in this function.
	// Inspect the first element to see if we have 12 parameters => proper motion fitting happens.
	double omega, e, tau, T, Omega, inc, alpha;
	double K, gamma;
	double x0, y0, mu_x, mu_y, pi;
	bool fit_pm = false;

	// Iterate over the possible modes stored in params and write out a unique residual data file for each mode.
	for(int mode = 0; mode < params.size(); mode++)
	{
		// Set the output filename, extract the parameters
		stringstream ast_output;
		ast_output << output_basename << "_"  << mode + 1 << "_res_ast";

		stringstream rv_output;
		rv_output << output_basename << "_"  << mode + 1 << "_res_rv";

		// Extract the parameters, then generate the data.
		extract_params_both(params[mode], fit_pm, omega, e, tau, T, K, gamma, Omega, inc, alpha, x0, y0, mu_x, mu_y, pi);

		residuals_ast(ast_output.str(), fit_pm, omega, e, tau, T, Omega, inc, alpha, x0, y0, mu_x, mu_y, pi, ast_data);
		residuals_rv(rv_output.str(), omega, e, tau, T, K, gamma, rv_data);
	}
}

void residuals_rv(string output_basename, vector<vector<double> > & params, vector<vector<double> > & data)
{
	// With the astrometry data we either fit 7, 8, 12, or 13 parameters in fitrv, but the noise parameter
	// doesn't change anything in this function.
	// Inspect the first element to see if we have 12 parameters => proper motion fitting happens.
	double omega, e, tau, T, K, gamma;

	// Iterate over the possible modes stored in params and write out a unique residual data file for each mode.
	for(int mode = 0; mode < params.size(); mode++)
	{
		// Set the output filename, extract the parameters
		stringstream output_name;
		output_name << output_basename << "_" << mode + 1 << "_res_rv";

		extract_params_rv(params[mode], omega, e, tau, T, K, gamma);

		residuals_rv(output_name.str(), omega, e, tau, T, K, gamma, data);
	}
}

void residuals_rv(string output_name,
		double omega, double e, double tau, double T,
		double K, double gamma,
		vector<vector<double> > & data)
{

	double n, M, E, cos_E, sin_E, beta;
	double sin_omega, cos_omega;
	double l1, m1, l2, m2, n1, n2;
	double t, rv, rvi;

	ofstream output;
	output.precision(10);
	output.setf(ios::fixed,ios::floatfield);
	output.open(output_name.c_str());

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

		rv = gamma + K * beta / (1 - e * cos_E) * (beta * cos_E * cos_omega - sin_E * sin_omega);

		output << t << " " << rv - rvi << endl;
	}

	output.close();

}

