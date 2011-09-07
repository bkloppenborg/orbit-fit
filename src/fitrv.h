/*
 * fitrv.h
 *
 *  Created on: Aug 15, 2011
 *      Author: bkloppenborg
 *
 *  Fit radial velocity data using a nested sampling Bayesian routine.
 */

#ifndef MAIN_H_
#define MAIN_H_

// Prints out help describing the options on the command line
void print_help();

// The main routine.  Basically just used to parse out some parameters before handing
// things off to other functions.
int main(int argc, char *argv[]);

#endif /* MAIN_H_ */
