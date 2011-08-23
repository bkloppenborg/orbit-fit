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
#define MAS_TO_RAD PI / 648000000
#define RAD_TO_MAS 648000000 / PI

#define TWO_PI 2*PI

#define DAY_IN_SEC 24 * 60 * 60
#define SEC_IN_DAY 1.0 / DAY_IN_SEC

double DegToRad(double value);

double MasToRad(double value);

double RadToDeg(double value);

double RadToMas(double value);

#endif //COMMON_H
