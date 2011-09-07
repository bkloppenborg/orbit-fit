/*
 * multinest_inf.h
 *
 *  Created on: Aug 22, 2011
 *      Author: bkloppenborg
 *
 *  Header for a common interface to multinest.
 */

#ifndef MULTINEST_INF_H_
#define MULTINEST_INF_H_

#include "multinest.h"

// A wrapper to kick off multinest.
void run_multinest(int mmodal, int ceff, int nlive, double tol,
    double efr, int ndims, int nPar, int nClsPar,
    int maxModes, int updInt, double Ztol, char root[],
    int seed, int *pWrap, int fb, int resume,
    int outfile, int initMPI, double logZero,
    void (*LogLike)(double *, int *, int *, double *),
    void (*dumper)(int *, int *,  int *, double **, double **, double *, double *, double *, double *),
    int context);

#endif /* MULTINEST_INF_H_ */
