/*
 * main.h
 *
 *  Created on: Aug 22, 2011
 *      Author: bkloppenborg
 *
 *  Fit astrometric data using a nested sampling Bayesian algorithm.
 */

#ifndef FITAST_COMMON_H
#define FITAST_COMMON_H

#include <vector>
#include <string>

#include "ReadTextFile.h"
#include "orbit.h"
#include "common.h"

using namespace std;

namespace fitast
{
	extern double s_min;
	extern double s_max;

	extern bool read_no_error;
	extern double default_error;

	void log_likelihood(double *Cube, int *ndim, int *npars, double *lnew);

	void dumper(int *nSamples, int *nlive, int *nPar, double **physLive, double **posterior, double *paramConstr, double *maxLogLike, double *logZ, double *logZerr);

	void compute_scales();
	void compute_partial_priors();

	void run_fit();

	void ParseProgOptions(int argc, char *argv[], bool & param_error);

	void print_param_limits();

}

#endif /* FITAST_COMMON_H */
