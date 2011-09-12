/*
 * common.cpp
 *
 *  Created on: Aug 15, 2011
 *      Author: bkloppenborg
 *
 *  Common routines used among the subprograms
 */

#include <cstring>
#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <vector>

#include "common.h"
#include "constants.h"

using namespace std;

vector<string> param_names;

double e_min;
double e_max;
double omega_min;
double omega_max;
double T_min;
double T_max;
double tau_min;
double tau_max;

double scale_e;
double scale_omega;
double scale_T;
double scale_tau;

double prior_e;
double prior_omega;
double prior_T;
double prior_tau;

double DegToRad(double value)
{
	return value * DEG_TO_RAD;
}

double MasToRad(double value)
{
	return value * MAS_TO_RAD;
}

double RadToDeg(double value)
{
	return value * RAD_TO_DEG;
}

double RadToMas(double value)
{
	return value * RAD_TO_MAS;
}

/* A common function for parsing command line options to the fitting routines.
*  Definition of input parameters:
*      Omega:  Position Angle (in radians)
*      inc:    Orbital Inclination (in radians) (I didn't use i here because i is often used for indexing.
*      omega:  Argument/longitude of periestron (in radians)
*      asini:  Semi-major axis, typically in linear units.
*      alpha:  Semi-major axis, angular units
*      e:      Orbital eccentricity
*      tau:    Time of periestron passage (units of time)
*      T:      Orbital Period (same units as tau)
*/
void ParseCommonParams(int argc, char *argv[], bool & param_error)
{
	// First initialize to broad defaults:
    omega_min = 0;
    omega_max = TWO_PI;
    e_min = 0;
    e_max = 1;
    tau_min = 2.3E6;
    tau_max = 2.5E6;
    T_min = 1;
    T_max = 1E4;

    // Parse common options to limit these values further.
	for (int i = 1; i < argc; i++)
	{

		if(strcmp(argv[i], "-omega_min") == 0)
		{
			omega_min = atof(argv[i + 1]) * DEG_TO_RAD;
		}

		if(strcmp(argv[i], "-omega_max") == 0)
		{
			omega_max = atof(argv[i + 1]) * DEG_TO_RAD;
		}

		if(strcmp(argv[i], "-e_min") == 0)
		{
			e_min = atof(argv[i + 1]);
		}

		if(strcmp(argv[i], "-e_max") == 0)
		{
			e_max = atof(argv[i + 1]);
		}

		if(strcmp(argv[i], "-tau_min") == 0)
		{
			tau_min = atof(argv[i + 1]);
		}

		if(strcmp(argv[i], "-tau_max") == 0)
		{
			tau_max = atof(argv[i + 1]);
		}

		if(strcmp(argv[i], "-T_min") == 0)
		{
			T_min = atof(argv[i + 1]);
		}

		if(strcmp(argv[i], "-T_max") == 0)
		{
			T_max = atof(argv[i + 1]);
		}

    }

	// Lastly check that *_min < *_max, print out an error message and quit
	if(omega_min > omega_max)
	{
		param_error = true;
		cout << "Error: omega_min exceeds omega_max, exiting.";
	}

	if(e_min > e_max)
	{
		param_error = true;
		cout << "Error: e_min exceeds e_max, exiting.";
	}

	if(tau_min > tau_max)
	{
		param_error = true;
		cout << "Error: tau_min exceeds tau_max, exiting.";
	}

	if(T_min > T_max)
	{
		param_error = true;
		cout << "Error: T_min exceeds T_max, exiting.";
	}

}

void print_common_param_limits()
{
	printf("omega (deg) : %1.4e %1.4e \n", omega_min * RAD_TO_DEG, omega_max * RAD_TO_DEG);
	printf("e           : %1.4e %1.4e \n", e_min, e_max);
	printf("T (time)    : %1.4e %1.4e \n", T_min, T_max);
	printf("tau (time)  : %1.4e %1.4e \n", tau_min, tau_max);
}
