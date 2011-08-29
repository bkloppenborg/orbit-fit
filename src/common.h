/*
 * common.h
 *
 *  Created on: Aug 15, 2011
 *      Author: bkloppenborg
 *
 *  Common routines used among the subprograms
 */

#ifndef COMMON_H
#define COMMON_H

#ifdef M_PI
#define PI M_PI
#else
#define PI 3.1415926535897932384626433832795028841968
#endif

#define DEG_TO_RAD PI / 180.0
#define RAD_TO_DEG 180.0 / PI
#define MAS_TO_RAD 4.8481368E-9
#define RAD_TO_MAS 1.0 / (4.8481368E-9)

#define TWO_PI 2*PI

#define DAY_TO_SEC (24 * 60 * 60)
#define SEC_TO_DAY 1.0 / (24*60*60)
#define YEAR_TO_SEC 31556926.
#define SEC_TO_YEAR 1.0 / (YEAR_TO_SEC)

double DegToRad(double value);

double MasToRad(double value);

double RadToDeg(double value);

double RadToMas(double value);

void ParseCommandLine(int argc, char *argv[],
		double & Omega_min, double & Omega_max, double & inc_min, double & inc_max,
		double & omega_min, double & omega_max, double & asini_min, double & asini_max,
		double & alpha_min, double & alpha_max, double & e_min, double & e_max,
		double & tau_min, double & tau_max, double & T_min, double & T_max,
		bool & param_error);

#endif //COMMON_H
