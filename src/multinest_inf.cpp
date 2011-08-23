/*
 * multinest_inf.cpp
 *
 *  Created on: Aug 22, 2011
 *      Author: bkloppenborg
 */

#include "multinest_inf.h"
#include <cstring>

using namespace std;

void dumper(int *nSamples, int *nlive, int *nPar, double **physLive,
    double **posterior, double *paramConstr, double *maxLogLike, double *logZ)
{
    // Do nothing.
}

// A wrapper to kick off multinest.
void run_multinest(int mmodal, int ceff, int nlive, double tol,
    double efr, int ndims, int nPar, int nClsPar,
    int maxModes, int updInt, double Ztol, char root[],
    int seed, int *pWrap, int fb, int resume,
    int outfile, int initMPI, double logZero,
    void (*LogLike)(double *, int *, int *, double *),
    void (*dumper)(int *, int *, int *, double **, double **, double *, double *, double *),
    int context)
{
    // Clear out the remaining characters in the string:
    int i;
	for (i = strlen(root); i < 100; i++)
	    root[i] = ' ';

    // Run the nested sampling algorithm
    NESTRUN(&mmodal, &ceff, &nlive, &tol,
        &efr, &ndims, &nPar, &nClsPar,
        &maxModes, &updInt, &Ztol, root,
        &seed, pWrap, &fb, &resume,
        &outfile, &initMPI, &logZero,
        LogLike,
        dumper,
        &context);
}
