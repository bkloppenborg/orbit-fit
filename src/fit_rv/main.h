/*
 * main.h
 *
 *  Created on: Aug 15, 2011
 *      Author: bkloppenborg
 */

#ifndef MAIN_H_
#define MAIN_H_

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

#include "ReadTextFile.h"
#include "orbit.h"
#include "common.h"
#include "multinest.h"

using namespace std;

// Prints out help describing the options on the command line
void print_help();

void read_data(string filename, string comment_chars, double defaut_error, vector< vector<int> > split_info, vector< vector<double> > & data);

void log_likelihood(double *Cube, int *ndim, int *npars, double *lnew);

void run_fit(vector< vector<double> > & data);

void dumper(int *nSamples, int *nlive, int *nPar, double **physLive,
    double **posterior, double *paramConstr, double *maxLogLike, double *logZ);

// A wrapper to kick off multinest.
void run_multinest(int mmodal, int ceff, int nlive, double tol,
    double efr, int ndims, int nPar, int nClsPar,
    int maxModes, int updInt, double Ztol, char root[],
    int seed, int *pWrap, int fb, int resume,
    int outfile, int initMPI, double logZero,
    void (*LogLike)(double *, int *, int *, double *),
    void (*dumper)(int *, int *, int *, double **, double **, double *, double *, double *),
    int context);

// The main routine.  Basically just used to parse out some parameters before handing
// things off to other functions.
int main(int argc, char *argv[]);

#endif /* MAIN_H_ */
