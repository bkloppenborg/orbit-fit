/*
 * main.h
 *
 *  Created on: Aug 15, 2011
 *      Author: bkloppenborg
 *
 *  Fit radial velocity data using a nested sampling Bayesian routine.
 */

#ifndef FITRV_COMMON_H_
#define FITRV_COMMON_H_

#include <string>
#include <vector>

using namespace std;

namespace fitrv
{

	void read_data(string filename, string comment_chars, vector< vector<int> > split_info);

	void dumper(int *nSamples, int *nlive, int *nPar, double **physLive, double **posterior, double *paramConstr, double *maxLogLike, double *logZ, double *logZerr);

	void log_likelihood(double *Cube, int *ndim, int *npars, double *lnew);

	void compute_scales();
	void compute_partial_priors();

	void run_fit();

	void ParseProgOptions(int argc, char *argv[], bool & param_error);

	void print_param_limits();

}

#endif /* FITRV_COMMON_H_ */
