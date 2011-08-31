/*
 * main.h
 *
 *  Created on: Aug 15, 2011
 *      Author: bkloppenborg
 *
 *  Fit radial velocity data using a nested sampling Bayesian routine.
 */

#ifndef MAIN_H_
#define MAIN_H_

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

#include "ReadTextFile.h"
#include "orbit.h"
#include "common.h"
#include "multinest.h"

using namespace std;

// Prints out help describing the options on the command line
void print_help();

void read_data(string filename, string comment_chars, double defaut_error, vector< vector<int> > split_info, vector< vector<double> > & data);

void dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double *paramConstr, double &maxLogLike, double &logZ, double &logZerr);

void log_likelihood(double *Cube, int *ndim, int *npars, double *lnew);

void run_fit(vector< vector<double> > & data);

void ParseProgOptions(int argc, char *argv[], bool & param_error);

// The main routine.  Basically just used to parse out some parameters before handing
// things off to other functions.
int main(int argc, char *argv[]);

#endif /* MAIN_H_ */
