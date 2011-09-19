/*
 * fitrv.cpp
 *
 *  Created on: Aug 14, 2011
 *      Author: bkloppenborg
 *
 *  Fit radial velocity data using a nested sampling Bayesian routine.
 */

#include "fitrv_common.h"

#include <string>
#include <iostream>
#include <cstdio>

#include <common.h>
#include "read_data.h"

using namespace std;
using namespace fitrv;

extern vector<string> param_names;

extern vector<vector<double> > rv_data;

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
    ParseCommonParams(argc, argv, param_error);
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

    read_data_rv(input_rv, comment_chars, split_info, rv_data, read_no_error, default_error);

    if(rv_data.size() == 0)
    {
    	printf("Data file is empty!  Exiting.\n");
    	return 0;
    }

    // Push the parameter names onto the name vector:
	param_names.push_back("omega ");
	param_names.push_back("e     ");
	param_names.push_back("T     ");
	param_names.push_back("tau   ");
	param_names.push_back("K     ");
	param_names.push_back("gamma ");
	param_names.push_back("s     ");

    run_fit();

	return 0;
}
