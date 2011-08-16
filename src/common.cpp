/*
 * common.cpp
 *
 *  Created on: Aug 15, 2011
 *      Author: bkloppenborg
 *
 *  Common routines used among the subprograms
 */

#include "common.h"

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
