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

double DegToRad(double value);

double MasToRad(double value);

double RadToDeg(double value);

double RadToMas(double value);

void ParseCommonParams(int argc, char *argv[], bool & param_error);

void print_common_param_limits();

#endif //COMMON_H
