// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_GEOM_TRAJ_CIRCULAR_H
#define XCOLL_GEOM_TRAJ_CIRCULAR_H

#include <stdio.h>
#include <math.h>


// /*gpufun*/
// void CircularTrajectory_set_params(CircularTrajectory traj, double sR, double xR,
//                                         LocalParticle part){
//     CircularTrajectory_set_sR(traj, sR);
//     CircularTrajectory_set_xR(traj, xR);
//     double s0 = LocalParticle_get_s(part);
//     double x0 = LocalParticle_get_x(part);
//     double R = sqrt((s0-sR)*(s0-sR) + (x0-xR)*(x0-xR));
//     CircularTrajectory_set_R(traj, R);
//     CircularTrajectory_set_sin_tI(traj, (x0-xR)/R);
//     CircularTrajectory_set_cos_tI(traj, (s0-sR)/R);
//     CircularTrajectory_set_tan_tI(traj, (x0-xR)/(s0-sR));
// }


// TODO: maybe for this trajectory it is faster to use tI directly instead of sin_tI and cos_tI

/*gpufun*/
void CircularTrajectory_set_params(CircularTrajectory traj, double sR, double xR,
                                     double s0, double x0){
    CircularTrajectory_set_sR(traj, sR);
    CircularTrajectory_set_xR(traj, xR);
    double R = sqrt((s0-sR)*(s0-sR) + (x0-xR)*(x0-xR));
    CircularTrajectory_set_R(traj, R);
    CircularTrajectory_set_sin_tI(traj, (x0-xR)/R);
    CircularTrajectory_set_cos_tI(traj, (s0-sR)/R);
    CircularTrajectory_set_tan_tI(traj, (x0-xR)/(s0-sR));
}

/*gpufun*/
double CircularTrajectory_func_s(CircularTrajectory traj, double l){
    double R = CircularTrajectory_get_R(traj);
    double sR = CircularTrajectory_get_sR(traj);
    double sin_tI = CircularTrajectory_get_sin_tI(traj);
    double cos_tI = CircularTrajectory_get_cos_tI(traj);
    return sR + R*cos(l)*cos_tI - R*sin(l)*sin_tI; // s(𝜆) = sR + R cos(𝜆 + 𝜃I)
}

/*gpufun*/
double CircularTrajectory_func_x(CircularTrajectory traj, double l){
    double R = CircularTrajectory_get_R(traj);
    double xR = CircularTrajectory_get_xR(traj);
    double sin_tI = CircularTrajectory_get_sin_tI(traj);
    double cos_tI = CircularTrajectory_get_cos_tI(traj);
    return xR + R*sin(l)*cos_tI + R*cos(l)*sin_tI; // s(𝜆) = sR + R sin(𝜆 + 𝜃I)
}

/*gpufun*/
double CircularTrajectory_func_xp(CircularTrajectory traj, double l){
    double tan_tI = CircularTrajectory_get_tan_tI(traj);
    return (tan_tI + tan(l)) / (1. - tan_tI*tan(l)); // 𝜃(𝜆) = 𝜃I + 𝜆 + chan. effects
}

/*gpufun*/
double CircularTrajectory_deriv_s(CircularTrajectory traj, double l){
    double R = CircularTrajectory_get_R(traj);
    double sin_tI = CircularTrajectory_get_sin_tI(traj);
    double cos_tI = CircularTrajectory_get_cos_tI(traj);
    return -R*sin(l)*cos_tI - R*cos(l)*sin_tI; // s(𝜆) = sR + R cos(𝜆 + 𝜃I)
}

/*gpufun*/
double CircularTrajectory_deriv_x(CircularTrajectory traj, double l){
    double R = CircularTrajectory_get_R(traj);
    double sin_tI = CircularTrajectory_get_sin_tI(traj);
    double cos_tI = CircularTrajectory_get_cos_tI(traj);
    return R*cos(l)*cos_tI - R*sin(l)*sin_tI; // s(𝜆) = sR + R sin(𝜆 + 𝜃I)
}

/*gpufun*/
void CircularTrajectory_bounding_box_s(CircularTrajectory traj, double l1, double l2, double extrema[2]){
    double s1 = CircularTrajectory_func_s(traj, l1);
    double s1_l0 = CircularTrajectory_func_s(traj, 0.0); 
    double s2 = CircularTrajectory_func_s(traj, l2);
    double sR = CircularTrajectory_get_sR(traj);
    double R  = CircularTrajectory_get_R(traj);
    extrema[0] = MIN(s1, s2);
    double theta = atan(CircularTrajectory_get_tan_tI(traj));
    if ((s1_l0 < sR)) {
        if (theta < 0){
            theta = M_PI + theta;
        } else {
            theta = -(M_PI/2 + (M_PI/2 - theta));
        }
    }
    double l1_rescaled = theta + l1;
    double l2_rescaled = l1_rescaled + fabs(l1-l2);
    if ((l1_rescaled <= M_PI && M_PI <= l2_rescaled) || (l1_rescaled <= 3*M_PI && 3*M_PI <= l2_rescaled)){
        extrema[0] = sR - R;
    } else {
        extrema[0] = MIN(s1, s2);
    }
    if ((l1_rescaled <= 0. && 0. <= l2_rescaled) || (l1_rescaled <= 2*M_PI && 2*M_PI <= l2_rescaled)){
        extrema[1] = sR + R;
    } else {
        extrema[1] = MAX(s1, s2);
    }
}

/*gpufun*/
void CircularTrajectory_bounding_box_x(CircularTrajectory traj, double l1, double l2, double extrema[2]){
    double x1 = CircularTrajectory_func_x(traj, l1);
    double x2 = CircularTrajectory_func_x(traj, l2);
    double R  = CircularTrajectory_get_R(traj);
    double xR = CircularTrajectory_get_xR(traj);
    double theta = atan(CircularTrajectory_get_tan_tI(traj));
    if ((CircularTrajectory_func_s(traj, 0.0) < CircularTrajectory_get_sR(traj))){
        if (theta < 0){
            theta = M_PI + theta;
        } else {
            theta = -(M_PI/2 + (M_PI/2 - theta));;
        }
    }
    double l1_rescaled = theta + l1;
    double l2_rescaled = l1_rescaled + fabs(l1-l2);
    extrema[0] = MIN(x1, x2);
    extrema[1] = MAX(x1, x2);
    if ((l1_rescaled <= 3*M_PI/2. && 3*M_PI/2. <= (l2_rescaled)) || (l1_rescaled <= -M_PI/2. && -M_PI/2. <= (l2_rescaled))){
        extrema[0] = xR - R;
    } 
    if ((l1_rescaled <= M_PI/2. && M_PI/2. <= (l2_rescaled)) || (l1_rescaled <= 5*M_PI/2. && 5*M_PI/2. <= (l2_rescaled))){
        extrema[1] = xR + R;
    }
}

#endif /* XCOLL_GEOM_TRAJ_CIRCULAR_H */