#include <math.h>
#include "orbit.h"
using namespace std;

/*
 * orbit.c
 *
 *  Created on: Aug 14, 2011
 *      Author: bkloppenborg
 *
 *  Implements basic routines for computing orbital solutions, given data and some
 *  fitting method (least squares, Bayesian, etc.)
 *   *
 *  All angles are in RADIANS!
 *
 *  Definition of input parameters:
 *      Omega:  Position Angle (in radians)
 *      inc:    Orbital Inclination (in radians) (I didn't use i here because i is often used for indexing.
 *      omega:  Argument/longitude of periestron (in radians)
 *      a:      Semi-major axis, typically in linear units.
 *      alpha:  Semi-major axis, angular units
 *      e:      Orbital eccentricity
 *      tau:    Time of periestron passage (units of time)
 *      T:      Orbital Period (same units as tau)
 *      t:		Time of observation (same units as tau)
 *
 *  Computed parameters:
 *      E:      Eccentric Anomaly (from Kepler's equation)
 *      M:      Mean anomaly
 *      n:      Mean angular velocity (2 pi / T)
 *
 *  Default ordering of parameters:
 *  (double Omega, double inc, double omega, a, e, tau, T)
 */

// TODO: It would be better to look this up in a library, but the precision here should be good enough for now.
#ifdef M_PI
#define PI M_PI
#else
#define PI 3.1415926535897932384626433832795028841968
#endif

// Computes the mean angular velocity for the given period, T
double ComputeN(double T)
{
	return 2 * PI / T;
}

// Computes the mean anomaly given tau, n, and the time of observation, t.
double ComputeM(double tau, double n, double t)
{
    return n*(t - tau);
}

double ComputeE(double M, double e)
{
    /*
    Solve Kepler's equation, M = E - e sin(E), using a Newton-Raphson method:

    TODO: This function should be rewritten to use Laguerre polynomials:
    http://www.springerlink.com/content/p122000960815647/
    as it has been shown to converge for all e, and converge faster than Newton-Raphson

    TODO: It would be a good idea to use memoization here, will speed up the minimization significantly.

    TODO: Tweak this value, along with the tolerance
    */

    int i = 0;

    //printf("M %f \n", M);
    // Initial guess (from Smith (1979))
    double E = M + e * sin(M) / (1 - sin(M + e) + sin(M));
    double E_old = 0;
    //printf("E0 %f \n", E);

    // Below we are using the Newton-Raphson solution.  The second term is the derivative
    // of Kepler's equation.
    if(e > 0.5)
    {
		do
		{
			i++;
			E_old = E;
			E = E_old + (M + e * sin(E_old) - E_old) / (1 - e * cos(E_old));
			//printf("E %f \n", E);
		}
		while ((i < MAX_ITERATIONS) && ((E - E_old) > THRESH));
    }
    else
    {
		do
		{
			i++;
			E_old = E;
			E = M + e * sin(E);
		} while ((i < MAX_ITERATIONS) && ((E - E_old) > THRESH));
    }

//    if(i == (MAX_ITERATIONS - 1))
//        printf("WARNING: Failed to Solve Kepler's Equation!");

    //printf("Ef %f \n", E);

    return E;
}

// Computes the positions (x, y, z) and velocities (v_x, v_y, v_z) at time t
// subject to the specified orbital parameters.
void GetAll(double Omega, double inc, double omega, double a, double e, double tau, double T, double t,
		double & x, double & y, double & z, double & v_x, double v_y, double & v_z)
{
	// Local variables:
    double l1, l2, m1, m2, n1, n2;

	// Pre-compute a few values
    double n = ComputeN(T);
    double M = ComputeM(tau, n, t);
    double E = ComputeE(M, e);

    double cos_E = cos(E);
    double sin_E = sin(E);

    // Now compute the orbital coefficients
    Compute_Coefficients(Omega, inc, omega, l1, m1, n1, l2, m2, n2);

    // For bound, elliptical orbits we have
    double beta = sqrt(1 - e*e);
    double eta = 1 - e * cos_E;

    // Compute the positions and velocities:
    Compute_xyz(a, beta, e, l1, l2, m1, m2, n1, n2, cos_E, sin_E, x, y, z);
    Compute_v_xyz(a, beta, e, l1, l2, m1, m2, n1, n2, cos_E, sin_E, n, eta, v_x, v_y, v_z);
}

void GetVelocities(double Omega, double inc, double omega, double a, double e, double tau, double T, double t,
	double & v_x, double & v_y, double & v_z)
{
	// Local variables
	double l1, l2, m1, m2, n1, n2;
    double n = ComputeN(T);
    double M = ComputeM(tau, n, t);
    double E = ComputeE(M, e);

    double cos_E = cos(E);
    double sin_E = sin(E);

    Compute_Coefficients(Omega, inc, omega, l1, m1, n1, l2, m2, n2);

    // For bound, elliptical orbits we have
    double beta = sqrt(1 - e*e);
    double eta = 1 - e * cos_E;

    // Compute the positions and velocities:
    Compute_v_xyz(a, beta, e, l1, l2, m1, m2, n1, n2, cos_E, sin_E, n, eta, v_x, v_y, v_z);
}

// Computes the radius vector and cross-checks it.
double ComputeR(double a, double e, double E)
{
    // First calculate r the easy way:
    double r1 = a * (1 - e * cos(E));

	// Now the slightly more involved method:
    double f = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2));
    double r2 = a * (1 - e*e) / (1 + e * cos(f));

//    if(fabs(r1 - r2) > THRESH)
//        printf("WARNING: RADII DO NOT AGREE!");

    return r2;
}


// Computes the coefficients, L1, L2, M1, M2, N1, N2 for the orbital equations
void Compute_Coefficients(double Omega, double inc, double omega,
		double & L1, double & M1, double & N1, double & L2, double & M2, double & N2)
{
    double c_Omega = cos(Omega);
    double s_Omega = sin(Omega);
    double c_inc = cos(inc);
    double s_inc = sin(inc);
    double c_omega = cos(omega);
    double s_omega = sin(omega);

    // Note, if an error is found here, be sure to update the GetRV and GetXY functions.
    L1 = c_Omega * c_omega - s_Omega * s_omega * c_inc;
    M1 = s_Omega * c_omega + c_Omega * s_omega * c_inc;
    N1 = s_omega * s_inc;
    L2 = -c_Omega * s_omega - s_Omega * c_omega * c_inc;
    M2 = -s_Omega * s_omega + c_Omega * c_omega * c_inc;
    N2 = c_omega * s_inc;
}

// Computes the (x, y, z) position of the orbit.
// Uses the equations defined in "Orbital Motions" by A. E. Roy 2005 pg. 93
void Compute_xyz(double a, double beta, double e,
		double l1, double l2, double m1, double m2, double n1, double n2,
		double cos_E, double sin_E,
		double & x, double & y, double & z)
{
    x = a * (l1 * cos_E + beta * l2 * sin_E - e * l1);
    y = a * (m1 * cos_E + beta * m2 * sin_E - e * m1);
    z = a * (n1 * cos_E + beta * n2 * sin_E - e * n1);
}

// Computes the (v_x, v_y, v_z) velocities of the orbit.  Equations defined in
// "Orbital Motions" by A. E. Roy 2005 pg. 93
void Compute_v_xyz(double a, double beta, double e,
		double l1, double l2, double m1, double m2, double n1, double n2,
		double cos_E, double sin_E, double n, double eta,
		double & v_x, double & v_y, double & v_z)
{
    v_x = n * a / eta * (beta * l2 * cos_E - l1 * sin_E);
    v_y = n * a / eta * (beta * m2 * cos_E - m1 * sin_E);
    v_z = n * a / eta * (beta * n2 * cos_E - n1 * sin_E);
}

// Computes the positions (x,y,z) for the given orbital elements at time t.
void GetPositions(double Omega, double inc, double omega, double a, double e, double tau, double T, double t,
	double & x, double & y, double & z)
{
	// Local variables:
    double l1, l2, m1, m2, n1, n2;

	// Pre-compute a few values
    double n = ComputeN(T);
    double M = ComputeM(tau, n, t);
    double E = ComputeE(M, e);

    double cos_E = cos(E);
    double sin_E = sin(E);

    // Now compute the orbital coefficients
    Compute_Coefficients(Omega, inc, omega, l1, m1, n1, l2, m2, n2);

    // For bound, elliptical orbits we have
    double beta = sqrt(1 - e*e);
    double eta = 1 - e * cos_E;

    // Compute the positions and velocities:
    Compute_xyz(a, beta, e, l1, l2, m1, m2, n1, n2, cos_E, sin_E, x, y, z);
}

// Computes only the radial velocity given the orbital parameters
// TODO: Add in a function that takes the more traditional K-based parameters too.
void GetRV(double inc, double omega, double a, double e, double tau, double T, double t,
	double & rv)
{
	// Local variables
    double n = ComputeN(T);
    double M = ComputeM(tau, n, t);
    double E = ComputeE(M, e);

    double cos_E = cos(E);
    double sin_E = sin(E);
    double s_inc = sin(inc);

    // For bound, elliptical orbits we have
    double beta = sqrt(1 - e*e);
    double eta = 1 - e * cos_E;

    double n1 = sin(omega) * s_inc;
    double n2 = cos(omega) * s_inc;

    // Compute the positions and velocities:
    rv = n * a / eta * (beta * n2 * cos_E - n1 * sin_E);
}

// Computes the radial velocity given the orbital parameters.  Uses asini instead of a and inc.
void GetRV(double omega, double asini, double e, double tau, double T, double t, double & rv)
{
	// Local variables
    double n = ComputeN(T);
    double M = ComputeM(tau, n, t);
    double E = ComputeE(M, e);

    double cos_E = cos(E);
    double sin_E = sin(E);

    // For bound, elliptical orbits we have
    double beta = sqrt(1 - e*e);
    double eta = 1 - e * cos_E;

    // Compute the positions and velocities:
    rv = n * asini / eta * (beta * cos(omega) * cos_E - sin(omega) * sin_E);
}
