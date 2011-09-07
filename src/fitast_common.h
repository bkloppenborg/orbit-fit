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
	void read_data(string filename, string comment_chars, vector< vector<int> > split_info);

	void log_likelihood(double *Cube, int *ndim, int *npars, double *lnew);

	void dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double *paramConstr, double &maxLogLike, double &logZ, double &logZerr);

	void run_fit();

	void ParseProgOptions(int argc, char *argv[], bool & param_error);

}

#endif /* FITAST_COMMON_H */
