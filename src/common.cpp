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
#include "common.h"

using namespace std;

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
int ParseCommandLine(int argc, char *argv[],
		double & Omega_min, double & Omega_max, double & inc_min, double & inc_max,
		double & omega_min, double & omega_max, double & asini_min, double & asini_max,
		double & alpha_min, double & alpha_max, double & e_min, double & e_max,
		double & tau_min, double & tau_max, double & T_min, double & T_max,
		bool & param_error)
{
	// First initialize to broad defaults:
	// NOTE: All angle units, except alpha, are converted to RADIANS before exiting this function
	// NOTE: alpha is kept in whatever native units are used for the input data
	// NOTE: All time units are converted to seconds
	Omega_min = 0;
	Omega_max = TWO_PI;
	inc_min = -90 * DEG_TO_RAD;
	inc_max = 90 * DEG_TO_RAD;
    omega_min = 0;
    omega_max = TWO_PI;
    alpha_min = 1;
    alpha_max = 100;
    asini_min = 1;
    asini_max = 1E20;
    e_min = 0;
    e_max = 1;
    tau_min = 2.5E3 * DAY_TO_SEC;
    tau_max = 2.5E6  * DAY_TO_SEC;
    T_min = 1 * DAY_TO_SEC;
    T_max = 1E5  * DAY_TO_SEC;

    // Parse common options to limit these values further.
	for (int i = 1; i < argc; i++)
	{
		if(strcmp(argv[i], "-Omega_min") == 0)
		{
			Omega_min = atof(argv[i + 1]) * DEG_TO_RAD;
		}

		if(strcmp(argv[i], "-Omega_max") == 0)
		{
			Omega_max = atof(argv[i + 1]) * DEG_TO_RAD;
		}

		if(strcmp(argv[i], "-inc_min") == 0)
		{
			inc_min = atof(argv[i + 1]) * DEG_TO_RAD;
		}

		if(strcmp(argv[i], "-inc_max") == 0)
		{
			inc_max = atof(argv[i + 1]) * DEG_TO_RAD;
		}

		if(strcmp(argv[i], "-omega_min") == 0)
		{
			omega_min = atof(argv[i + 1]) * DEG_TO_RAD;
		}

		if(strcmp(argv[i], "-omega_max") == 0)
		{
			omega_max = atof(argv[i + 1]) * DEG_TO_RAD;
		}

		if(strcmp(argv[i], "-asini_min") == 0)
		{
			asini_min = atof(argv[i + 1]);
		}

		if(strcmp(argv[i], "-asini_max") == 0)
		{
			asini_max = atof(argv[i + 1]);
		}

		if(strcmp(argv[i], "-alpha_min") == 0)
		{
			alpha_min = atof(argv[i + 1]);
		}

		if(strcmp(argv[i], "-alpha_max") == 0)
		{
			alpha_max = atof(argv[i + 1]);
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
			tau_min = atof(argv[i + 1]) * DAY_TO_SEC;
		}

		if(strcmp(argv[i], "-tau_max") == 0)
		{
			tau_max = atof(argv[i + 1]) * DAY_TO_SEC;
		}

		if(strcmp(argv[i], "-T_min") == 0)
		{
			T_min = atof(argv[i + 1]) * DAY_TO_SEC;
		}

		if(strcmp(argv[i], "-T_max") == 0)
		{
			T_max = atof(argv[i + 1]) * DAY_TO_SEC;
		}

		// Lastly check that *_min < *_max, print out an error message and quit
		if(Omega_min > Omega_max)
		{
			param_error = true;
			cout << "Error: Omega_min exceeds Omega_max.";
		}

		if(inc_min > inc_max)
		{
			param_error = true;
			cout << "Error: inc_min exceeds inc_max, exiting.";
		}

		if(omega_min > omega_max)
		{
			param_error = true;
			cout << "Error: omega_min exceeds omega_max, exiting.";
		}

		if(asini_min > asini_max)
		{
			param_error = true;
			cout << "Error: asini_min exceeds asini_max, exiting.";
		}

		if(alpha_min > alpha_max)
		{
			param_error = true;
			cout << "Error: alpha_min exceeds alpha_max, exiting.";
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
}
