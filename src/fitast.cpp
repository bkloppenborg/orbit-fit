/*
 * main.cpp
 *
 *  Created on: Aug 14, 2011
 *      Author: bkloppenborg
 *
 *  Fit astrometric data using a nested sampling Bayesian algorithm.
 */

#include <string>
#include <vector>
#include <iostream>

#include "fitast_common.h"
#include "read_data.h"

using namespace std;
using namespace fitast;

extern vector< vector<double> > ast_data;
extern bool fitast::read_no_error;
extern double fitast::default_error;


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
    ParseCommonParams(argc, argv, param_error);
    ParseProgOptions(argc, argv, param_error);

    if(param_error)
    	return 0;

    // Read in the file of RV data.
    // Each row should contain time, RV, [errors]
    const string comment_chars("\\#~$&£%");

    vector< vector<int> > split_info;
    read_data_ast(input_rv, comment_chars, split_info, ast_data, fitast::read_no_error, fitast::default_error);

    run_fit();

	return 0;
}
