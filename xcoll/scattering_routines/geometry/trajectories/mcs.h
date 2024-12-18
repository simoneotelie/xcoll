// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_COLL_GEOM_TRAJECTORIES_MCS_H
#define XCOLL_COLL_GEOM_TRAJECTORIES_MCS_H

#include "../c_init/simpson.h" // this is obv wrong 
#include <stdio.h>
#include <math.h>

/*gpufun*/
// int8_t DriftTrajectory_vlimit(double* restrict_s, double s0, double y0, double ym, double ymin, double ymax){
//     if (fabs(ym) < XC_EPSILON){
//         // Trajectory parallel to s axis
//         if (y0 < ymin || y0 > ymax){
//             return 0;  // Completely outside - no crossing possible
//         } else {
//             return -1; // Completely inside - no vertical check needed
//         }
//     } else {
//         restrict_s[0] = (ymin - y0)/ym + s0;
//         restrict_s[1] = (ymax - y0)/ym + s0;
//         SWAP(restrict_s, 0, 1);   // To make sure these are sorted
//         return 1;  // Default behavior: check overlap with horizontal crossings
//     }
// }

// MULTIPLE COULOMB SCATTERING TRAJECTORY ------------------------------------------------------------------

void generate_gaussian_random_numbers(double* z1, double* z2) {
    // get z1 and z2 from Box-Muller transform
    double v1, v2, r2;
    do {
        v1 = 2 * ((double)rand() / RAND_MAX) - 1;
        v2 = 2 * ((double)rand() / RAND_MAX) - 1;
        r2 = v1 * v1 + v2 * v2;
    } while (r2 >= 1 || r2 == 0);
    double a = sqrt(-2 * log(r2) / r2);
    *z1 = v1 * a;
    *z2 = v2 * a;
}

void A(double Xo, double p, double* Ax, double* Ay) {
    // z1, z2 are gaussian random numbers. Xo is the radiation length, p is the momentum
    double z1, z2;
    generate_gaussian_random_numbers(&z1, &z2);
    *Ax = (z1 / sqrt(12) + z2 / 2.0) * Xo * 13.6e-3 / p * 0.038;
    *Ay = (z1 / sqrt(12) + z2 / 2.0) * Xo * 13.6e-3 / p * 0.038;
}

// void B(double* B) {
//     *B = 1.0/0.038// + log(q**2/beta**2) // do we approx beta to one? 
// }

double MultipleCoulombTrajectory(double s, double const Xo, const double Ax, const double Ay, const double B){
    // MCS trajectory form PDG rewritted in terms of A, B and s/Xo. 
    return Ax * pow(sqrt(s/Xo),3.0) * (1.0/0.038 + log(s/Xo));
}


// CURVE LENGTH ---------------------------------------------------------------------------------------------
// Curve length is int_s1^s2 sqrt(1 + (dx/ds)^2 + (dy/ds)^2) ds

double integrand(double s, void *params) {
    // Extract Ax, Ay, Xo from the params
    double *p = (double *)params;
    const double Ax = p[0];
    const double Ay = p[1];
    const double Xo = p[2];
    double x;
    double y;

    x = Ax/Xo * (sqrt(s/Xo)*3.0/2.0*log(s/Xo)+1.0/0.038 + sqrt(s/Xo));
    y = Ay/Xo * (sqrt(s/Xo)*3.0/2.0*log(s/Xo)+1.0/0.038 + sqrt(s/Xo));
    return sqrt(1 + x*x + y*y);
}

/*gpufun*/
double MultipleCoulombTrajectory_length(double s0, double x0, double y0, const double s1, const double s2, 
                                        const double* Ax, const double* Ay, const double Xo){
    (void) s0;  // Avoid unused parameter warning
    (void) x0;  // Avoid unused parameter warning
    (void) y0;  // Avoid unused parameter warning
    double n = 1000; // number of subintervals, must be even
    double params[3] = {*Ax, *Ay, Xo};

    // compare with adaptive gauss/kronrod quadrature later! 
    double length = simpson(integrand, s1, s2, n, params); // shold avoid magic numbers
    return length;
}
