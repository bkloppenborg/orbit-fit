/*
 * orbit.h
 *
 *  Created on: Aug 14, 2011
 *      Author: bkloppenborg
 *
 *  Header file for orbit.c:
 *
 *  Implements basic routines for computing orbital solutions, given data and some
 *  fitting method (least squares, Bayesian, etc.)
 */

#ifndef ORBIT_H
#define ORBIT_H

// Define a few constants:
double const THRESH = 1E-5;
int const MAX_ITERATIONS = 50;

// Computes the mean angular velocity for the given period, T
double ComputeN(double T);

// Computes the mean anomaly given tau, n, and the time of observation, t.
double ComputeM(double tau, double n, double t);

double ComputeE(double M, double e);

// Computes the positions (x, y, z) and velocities (v_x, v_y, v_z) at time t
// subject to the specified orbital parameters.
void GetAll(double Omega, double inc, double omega, double a, double e, double tau, double T, double t,
		double & x, double & y, double & z, double & v_x, double v_y, double & v_z);

void GetVelocities(double Omega, double inc, double omega, double a, double e, double tau, double T, double t,
	double & v_x, double & v_y, double & v_z);

// Computes the radius vector and cross-checks it.
double ComputeR(double a, double e, double E);

// Computes the coefficients, L1, L2, M1, M2, N1, N2 for the orbital equations
void Compute_Coefficients(double Omega, double inc, double omega,
		double & L1, double & M1, double & N1, double & L2, double & M2, double & N2);

// Computes the (x, y, z) position of the orbit.
// Uses the equations defined in "Orbital Motions" by A. E. Roy 2005 pg. 93
void Compute_xyz(double a, double beta, double e,
		double l1, double l2, double m1, double m2, double n1, double n2,
		double cos_E, double sin_E,
		double & x, double & y, double & z);

// Computes the (v_x, v_y, v_z) velocities of the orbit.  Equations defined in
// "Orbital Motions" by A. E. Roy 2005 pg. 93
void Compute_v_xyz(double a, double beta, double e,
		double l1, double l2, double m1, double m2, double n1, double n2,
		double cos_E, double sin_E, double n, double eta,
		double & v_x, double & v_y, double & v_z);

// Computes the positions (x,y,z) for the given orbital elements at time t.
void GetPositions(double Omega, double inc, double omega, double a, double e, double tau, double T, double t,
	double & x, double & y, double & z);

// Computes only the radial velocity given the orbital parameters
// TODO: Add in a function that takes the more traditional K-based parameters too.
void GetRV(double inc, double omega, double a, double e, double tau, double T, double t,
	double & rv);

#endif // ORBIT_H
