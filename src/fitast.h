/*
 * fitast.h
 *
 *  Created on: Aug 22, 2011
 *      Author: bkloppenborg
 *
 *  Fit astrometric data using a nested sampling Bayesian algorithm.
 */

#ifndef FITAST_H_
#define FITAST_H_

#include "fitast_common.h"

using namespace std;

// Prints out help describing the options on the command line
void print_help();

// The main routine.  Basically just used to parse out some parameters before handing
// things off to other functions.
int main(int argc, char *argv[]);

#endif /* FITAST_H_ */
