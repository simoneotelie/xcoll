// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #
#include <math.h>
#ifndef XCOLL_GEOM__SEG_CIRCULAR_H
#define XCOLL_GEOM__SEG_CIRCULAR_H


/*gpufun*/
int8_t CircularSegment_bounded_below(CircularSegment seg){
    return 1;  // closed segment
}

/*gpufun*/
int8_t CircularSegment_bounded_above(CircularSegment seg){
    return 1;  // closed segment
}

/*gpufun*/
double CircularSegment_func_s(CircularSegment seg, double t){
    // t in [0,1] becomes tt in [theta1, theta2]
    double R = CircularSegment_get_R(seg);
    double sR = CircularSegment_get_sR(seg);
    double theta1 = CircularSegment_get__theta1(seg);
    double theta2 = CircularSegment_get__theta2(seg);
    return sR + R*cos(theta1 + t*(theta2 - theta1));
}

/*gpufun*/
double CircularSegment_func_x(CircularSegment seg, double t){
    double R = CircularSegment_get_R(seg);
    double xR = CircularSegment_get_xR(seg);
    double theta1 = CircularSegment_get__theta1(seg);
    double theta2 = CircularSegment_get__theta2(seg);
    return xR + R*sin(theta1 + t*(theta2 - theta1));
}

/*gpufun*/
double CircularSegment_deriv_s(CircularSegment seg, double t){
    double R = CircularSegment_get_R(seg);
    double theta1 = CircularSegment_get__theta1(seg);
    double theta2 = CircularSegment_get__theta2(seg);
    return -R*sin(theta1 + t*(theta2 - theta1))*(theta2 - theta1);
}

/*gpufun*/
double CircularSegment_deriv_x(CircularSegment seg, double t){
    double R = CircularSegment_get_R(seg);
    double theta1 = CircularSegment_get__theta1(seg);
    double theta2 = CircularSegment_get__theta2(seg);
    return R*cos(theta1 + t*(theta2 - theta1))*(theta2 - theta1);
}

/*gpufun*/
void CircularSegment_bounding_box(CircularSegment seg, double t1, double t2, BoundingBox* box){
    double s1 = CircularSegment_func_s(seg, t1);
    double s2 = CircularSegment_func_s(seg, t2);
    double x1 = CircularSegment_func_x(seg, t1);
    double x2 = CircularSegment_func_x(seg, t2);
    double sR = CircularSegment_get_sR(seg);
    double R  = CircularSegment_get_R(seg);
    double theta1 = CircularSegment_get__theta1(seg);
    double theta2 = CircularSegment_get__theta2(seg);
    double tt1 = theta1 + t1*(theta2 - theta1); // rescale
    double tt2 = theta1 + t2*(theta2 - theta1); // rescale
    double sin_t = (x2 - x1) / sqrt((x2 - x1)*(x2 - x1) + (s2 - s1)*(s2 - s1)); // box angle
    double cos_t = (s2 - s1) / sqrt((x2 - x1)*(x2 - x1) + (s2 - s1)*(s2 - s1));
    double sin_p, cos_p;
    if (sin_t < 0){   // if theta is larger than 180 degrees, theta = theta - 180
        sin_t = -sin_t;
        cos_t = -cos_t;
    }
    if (cos_t < 1){   // if theta is larger than 90 degrees, phi = theta + 90 
        sin_p = cos_t;
        cos_p = -sin_t;
    } else {          // if theta is between 0 and 90 degrees, phi = theta - 90
        sin_p = -cos_t;
        cos_p = sin_t;
    }
    if ( (tt2 - tt1) < 180){ // delta theta of box less than 180.
        box->l  = sqrt((s2 - s1)*(s2 - s1) + (x2 - x1)*(x2 - x1));   // length of the box
        box->w  = R - sqrt(R*R - box->l*box->l/4.);                  // width of the box, Sagitta
        // finding the first vertex. 
        if (x1 < x2){
            if (((t1 < M_PI && t2 > M_PI) && (cos_t < 0)) || ((-M_PI < t1 < 0.) && (-M_PI < t2 < 0.)) || ((t1 < 0 && t2 > 0) && (cos_t > 0))){                // t1 or t2 is lower vertex
                box->rC = sqrt( (s1+box->w*cos_p)*(s1+box->w*cos_p) +    // length of the position vector to the first vertex
                                (x1+box->w*sin_p)*(x1+box->w*sin_p) );
            } else {
                box->rC = sqrt(s1*s1 + x1*x1); // length of position vector to first vertex
            }
        } else {
            if (((t1 < M_PI && t2 > M_PI) && (cos_t < 0)) || ((-M_PI < t1 < 0.) && (-M_PI < t2 < 0.)) || ((t1 < 0 && t2 > 0) && (cos_t > 0))){                // t1 or t2 is lower vertex
                box->rC = sqrt( (s2+box->w*cos_p)*(s2+box->w*cos_p) +    // length of the position vector to the first vertex
                                (x2+box->w*sin_p)*(x2+box->w*sin_p) );
            } else {
                box->rC = sqrt(s2*s2 + x2*x2); // length of position vector to first vertex
            }
        }
    } else {
        double minimum, maximum;
        if (tt1 <= M_PI && M_PI <= tt2){
            minimum = sR - R;
        } else {
            minimum = MIN(s1, s2);
        }
        if ((tt1 <= 0. && 0. <= tt2) || (tt1 <= 2*M_PI && 2*M_PI <= tt2)){
            maximum = sR + R;
        } else {
            maximum = MAX(s1, s2);
        }
        box->rC = sqrt(

}

/*gpufun*/
void CircularSegment_bounding_box_x(CircularSegment seg, double t1, double t2, double extrema[2]){
    double x1 = CircularSegment_func_x(seg, t1);
    double x2 = CircularSegment_func_x(seg, t2);
    double R  = CircularSegment_get_R(seg);
    double xR = CircularSegment_get_xR(seg);
    double theta1 = CircularSegment_get__theta1(seg);
    double theta2 = CircularSegment_get__theta2(seg);
    double tt1 = theta1 + t1*(theta2 - theta1); // rescale
    double tt2 = theta1 + t2*(theta2 - theta1); // rescale
    if ((tt1 <= -M_PI/2. && -M_PI/2. <= tt2) || (tt1 <= 3*M_PI/2. && 3*M_PI/2. <= tt2)){
        extrema[0] = xR - R;
    } else {
        extrema[0] = MIN(x1, x2);
    }
    if (tt1 <= M_PI/2. && M_PI/2. <= tt2){
        extrema[1] = xR + R;
    } else {
        extrema[1] = MAX(x1, x2);
    }
}


}


// /*gpufun*/
// void CircularSegment_crossing_drift(CircularSegment seg, int8_t* n_hit, double* s, double s0, double x0, double xm){
//     // Get segment data
//     double R  = CircularSegment_get_R(seg);
//     double sC = CircularSegment_get_s(seg);
//     double xC = CircularSegment_get_x(seg);
//     double t1 = CircularSegment_get_t1(seg);
//     double t2 = CircularSegment_get_t2(seg);
//     // Move the angles to [-pi, pi]
//     int8_t reversed = 0, full_circle = 0;
//     if (fabs(fabs(t2 - t1) - 2*M_PI) < XC_EPSILON){
//         full_circle = 1;
//     }
//     while (t1 < -M_PI){
//         t1 += 2*M_PI;
//     }
//     while (t1 > M_PI){
//         t1 -= 2*M_PI;
//     }
//     while (t2 < -M_PI){
//         t2 += 2*M_PI;
//     }
//     while (t2 > M_PI){
//         t2 -= 2*M_PI;
//     }
//     if (t2 < t1){
//         reversed = 1;
//     }
//     // Calculate crossings
//     double a = 1 + xm*xm;
//     double bb = sC - xm*(x0 - xC - xm*s0); // This is -b/2 with b from the quadratic formula
//     double c = sC*sC + (x0 - xC - xm*s0)*(x0 - xC - xm*s0) - R*R;
//     double disc = bb*bb - a*c; // This is  2*discriminant**2
//     if (disc < 0){
//         // No crossing
//         return;
//     }
//     for (int8_t i = 0; i < 2; i++) {
//         double sgnD = i*2-1; // negative and positive solutions; if multiplicity 2, we add the same solution twice
//         double new_s = (bb + sgnD*sqrt(fabs(disc)))/a;
//         double new_x = x0 + (new_s - s0)*xm;
//         double t = atan2(new_x - xC, new_s - sC);
//         if (full_circle){
//             // Full circle, so always hit
//             s[*n_hit] = new_s;
//             (*n_hit)++;
//         } else if (reversed){
//             // t2 < t1, so we are looking at the inverted region of angles
//             if (t1 <= t || t <= t2){
//                 s[*n_hit] = new_s;
//                 (*n_hit)++;
//             }
//         } else {
//             if (t1 <= t && t <= t2){
//                 s[*n_hit] = new_s;
//                 (*n_hit)++;
//             }
//         }
//     }
// }
// /*gpufun*/
// void CircularSegment_crossing_mcs(CircularSegment seg, int8_t* n_hit, double* s, double x, const double* Ax, const double Xo){
//     return grid_search_and_newton();
// }
#endif /* XCOLL_GEOM__SEG_CIRCULAR_H */