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
void CircularTrajectory_init_bounding_box(CircularTrajectory traj, BoundingBox box, double l1, double l2){
    double theta = atan2(CircularTrajectory_get_sin_tI(traj), CircularTrajectory_get_cos_tI(traj));
    double ll1 = theta + l1;
    double ll2 = ll1 + (l2-l1);
    double s1 = CircularTrajectory_func_s(traj, l1);
    double x1 = CircularTrajectory_func_x(traj, l1);
    printf("\n\n\nNEW BOX\n");
    printf("s1 = %f, x1 = %f\n", s1, x1);
    printf("l1 = %f, l2 = %f\n", l1, l2);
    printf("theta = %f\n", theta);
    printf("ll1 = %f, ll2 = %f\n", ll1, ll2);
    double s2 = CircularTrajectory_func_s(traj, l2);
    double x2 = CircularTrajectory_func_x(traj, l2);
    double sR = CircularTrajectory_get_sR(traj);
    double xR = CircularTrajectory_get_xR(traj);
    double R  = CircularTrajectory_get_R(traj); 
    double dx = x2 - x1;
    double ds = s2 - s1;
    double chord_length = sqrt(dx*dx + ds*ds);
    double sin_chord = dx / chord_length;
    double cos_chord = ds / chord_length;
    double sin_t, cos_t;   // angle of the box wrt horizontal
    double min_x, min_s;
    double sin_rot, cos_rot;
    int8_t sign = 1;
    printf("cos_chord = %.10f, sin_chord = %f\n", cos_chord, sin_chord);
    if ((cos_chord > 1e-10) && (sin_chord > 1e-10)){    // if 0 < chord angle < 90 deg, then chord angle = box angle
        sin_t = sin_chord;
        cos_t = cos_chord;
        sin_rot = -cos_t;
        cos_rot = sin_t;
        sign = -1;
        printf("cos_t = %.10f, sin_t = %f\n", cos_t, sin_t);
        printf("cos_rot = %.10f, sin_rot = %f\n", cos_rot, sin_rot);
        printf("Je suis ici\n");
    } else {
        if (sin_chord < 1e-10) {   // if theta is larger than 180 degrees, theta = theta - 180
            sin_chord = -sin_chord;
            cos_chord = -cos_chord;
            printf("cos_t = %.10f, sin_t = %f\n", cos_chord, sin_chord);
            printf("Je suis ici 2\n");
        }
        if ((cos_chord > 1e-10) && (sin_chord > 1e-10)){    // if 0 < chord angle < 90 deg, then chord angle = box angle
            sin_t = sin_chord;
            cos_t = cos_chord;
            sin_rot = -cos_t;
            cos_rot = sin_t;
            sign = -1;
            printf("cos_t = %.10f, sin_t = %f\n", cos_t, sin_t);
            printf("cos_rot = %.10f, sin_rot = %f\n", cos_rot, sin_rot);
            printf("Je suis ici 3\n");
        } else {
            sin_t = - cos_chord;   // box angle is 90 degrees less than chord angle
            cos_t = sin_chord;
            sin_rot = sin_t;    
            cos_rot = cos_t;
            sign = -1;
            printf("cos_t = %.10f, sin_t = %f\n", cos_t, sin_t);
            printf("cos_rot = %.10f, sin_rot = %f\n", cos_rot, sin_rot);
            printf("Je suis ici 4\n");
        }
        printf("cos_chord = %.10f, sin_chord = %f\n", cos_chord, sin_chord);   
    }
    if ((ll2 - ll1) < M_PI){ // delta theta of box less than 180.
        BoundingBox_set_l(box, chord_length);   // length of the box
        BoundingBox_set_w(box, R - sqrt(R*R - chord_length*chord_length/4.));                  // width of the box, Sagitta
        double w = BoundingBox_get_w(box);
        // finding the first vertex. 
        if (x1 < x2){
            if (((l1 < M_PI && l2 > M_PI) && (cos_chord < 0)) || (( (-M_PI < l1) && (l1 < 0.) ) && ( (-M_PI < l2) && (l2 < 0.))) || ((l1 < 0 && l2 > 0) && (cos_chord > 0))){                // l1 or l2 is lower vertex
                BoundingBox_set_rC(box, sqrt( (s1+w*cos_t)*(s1+w*cos_t) +    // length of the position vector to the first vertex
                                              (x1+w*sin_t)*(x1+w*sin_t)) );
            } else {
                BoundingBox_set_rC(box, sqrt(s1*s1 + x1*x1)); // length of position vector to first vertex
            }
        } else {
            if (((l1 < M_PI && l2 > M_PI) && (cos_chord < 0)) || (( (-M_PI < l1) && (l1 < 0.) ) && ( (-M_PI < l2) && (l2 < 0.))) || ((l1 < 0 && l2 > 0) && (cos_chord > 0))){                // l1 or l2 is lower vertex
                BoundingBox_set_rC(box, sqrt( (s2+w*cos_t)*(s2+w*cos_t) +    // length of the position vector to the first vertex
                                              (x2+w*sin_t)*(x2+w*sin_t)) );
            } else {
                BoundingBox_set_rC(box, sqrt(s2*s2 + x2*x2)); // length of position vector to first vertex
            }
        }
    } else {
        double chord_side = 1;
        double l = (2*R);                               // L is always on the side of the chord
        double w = (R + sqrt(R*R - (chord_length*chord_length/4.))); // 2R - sagitta
        printf("2 r = %f\n", 2*R);
        printf("sagitta = %f\n", R - sqrt(R*R - chord_length*chord_length/4.));
        double chord_side_w1_x, chord_side_w1_s, chord_side_w2_x, chord_side_w2_s;
        double w3_x, w3_s, w4_x, w4_s;
        printf("w = %f\n", w);
        printf("l = %f\n", l);
        printf("chord length = %f\n", chord_length);
        printf("sin_chord = %f, cos_chord = %f\n", sin_chord, cos_chord);
        printf("sin_rot = %f, cos_rot = %f\n", sin_rot, cos_rot);
        printf("x1 = %f, x2 = %f\n", x1, x2);
        // determining the (s,x) of the vertices on the chord side of box
        if (x1 < x2){
            chord_side = -1;
            printf("!!!!!!!!!!!!! X1 < X2\n");
            printf("sign = %d\n", sign);
            printf("sin_t = %f, cos_t = %f\n", sin_t, cos_t);
            printf("s1 = %f, s2 = %f\n", s1, s2);
            printf("x1 = %f, x2 = %f\n", x1, x2);
            printf("s2 -2r - chord_length cos - = %f\n", s2 - (2*R - chord_length)/2.*cos_chord);
            printf("s2 + 2r - chord_length cos + = %f\n", s2 + (2*R - chord_length)/2.*cos_chord);
            printf("new point with cos = %f\n", (2*R - chord_length)/2.*cos_chord);
            printf("new point with sin = %f\n", (2*R - chord_length)/2.*sin_chord);
            printf("x2 - 2r - chord_length /2 sin -= %f\n", x2 - (2*R - chord_length)/2.*sin_chord);
            printf("x2 + 2r - chord_length /2 sin += %f\n", x2 + (2*R - chord_length)/2.*sin_chord);
            printf("s1 - 2r - chord_length /2 cos -= %f\n", s1 - (2*R - chord_length)/2.*cos_chord);
            printf("s1 + 2r - chord_length /2 cos += %f\n", s1 + (2*R - chord_length)/2.*cos_chord);
            printf("x1 - 2r - chord_length /2 sin -= %f\n", x1 - (2*R - chord_length)/2.*sin_chord);
            printf("x1 + 2r - chord_length /2 sin += %f\n", x1 + (2*R - chord_length)/2.*sin_chord);
            chord_side_w1_x = x1 - (2*R - chord_length)/2.*sin_chord;
            chord_side_w1_s = s1 - (2*R - chord_length)/2.*cos_chord;
            chord_side_w2_x = x2 + (2*R - chord_length)/2.*sin_chord;
            chord_side_w2_s = s2 + (2*R - chord_length)/2.*cos_chord;
            w3_x = chord_side_w1_x - sign*w*sin_rot; // sign because of the orientation of the box
            w3_s = chord_side_w1_s + w*cos_rot;
            w4_x = chord_side_w2_x - sign*w*sin_rot;
            w4_s = chord_side_w2_s + w*cos_rot;
            printf("w3: w1x - sign*w*sin_t = %f\n", chord_side_w1_x - sign*w*sin_rot);
            printf("w1x + w*sin_rot = %f\n", chord_side_w1_x + sign*w*sin_rot);
            printf("w2x - sign*w*sin_rot = %f\n", chord_side_w2_x - sign*w*sin_rot);
            printf("w2x + w*sin_rot = %f\n", chord_side_w2_x + sign*w*sin_rot);

            printf("cos_ti = %f, sin_ti = %f\n", cos_t, sin_t);
            printf("chord_side_w1_x = %f, chord_side_w1_s = %f\n", chord_side_w1_x, chord_side_w1_s);
        } else {
            chord_side = 1;
            printf("x1 = %f, x2 = %f\n", x1, x2);
            printf("s1 = %f, s2 = %f\n", s1, s2);
            printf("chord_length = %f\n", chord_length);
            printf("2r - chord_length cos - = %f\n", s2 - (2*R - chord_length)/2.*cos_chord);
            printf("2r - chord_length cos + = %f\n", s2 + (2*R - chord_length)/2.*cos_chord);
            printf("same without s2, cos = %f\n", (2*R - chord_length)/2.*cos_chord);
            printf("2r - chord_length/2 sin = %f\n", (2*R - chord_length)/2.*sin_chord);
            printf("2r - chord_length /2 sin -= %f\n", x2 - (2*R - chord_length)/2.*sin_chord);
            printf("2r - chord_length /2 sin += %f\n", x2 + (2*R - chord_length)/2.*sin_chord);
            printf("sint = %f, cost = %f\n", sin_t, cos_t);
            printf("L SIN = %f\n", (2*R - chord_length)/2.*sin_t);
            printf("L COS = %f\n", (2*R - chord_length)/2.*cos_t);
            chord_side_w2_x = x2 - (2*R - chord_length)/2.*sin_chord;
            chord_side_w2_s = s2 - (2*R - chord_length)/2.*cos_chord;
            chord_side_w1_x = x1 + (2*R - chord_length)/2.*sin_chord;
            chord_side_w1_s = s1 + (2*R - chord_length)/2.*cos_chord;
            w3_x = chord_side_w1_x + sign*w*sin_rot; 
            w3_s = chord_side_w1_s - w*cos_rot;
            w4_x = chord_side_w2_x + sign*w*sin_rot;
            w4_s = chord_side_w2_s - w*cos_rot;
            printf("cos_ti = %f, sin_ti = %f\n", cos_t, sin_t);
            printf("chord_side_w1_x = %f, chord_side_w1_s = %f\n", chord_side_w1_x, chord_side_w1_s);
        }
        printf("sign = %d\n", sign);
        printf("w3_x = %f, w3_s = %f\n", w3_x, w3_s);
        printf("w4_x = %f, w4_s = %f\n", w4_x, w4_s);
        printf("chord_side_w2_x = %f, chord_side_w2_s = %f\n", chord_side_w2_x, chord_side_w2_s);
        // Compare with the other three points
        min_x = chord_side_w1_x;
        min_s = chord_side_w1_s;
        // finding the first vertex. Also taking into account if the box is horiztonal
        printf("chordside = %f\n", chord_side);
        if ((chord_side_w2_x < min_x) || (chord_side_w2_x == min_x && chord_side_w2_s < min_s)) {
            if ((chord_side_w2_x < min_x) && (chord_side_w2_s > min_s)) {
                chord_side = -1;
            } else {
                chord_side = 1;
            }
            min_x = chord_side_w2_x; 
            min_s = chord_side_w2_s;
        }
        if (w3_x < min_x || (w3_x == min_x && w3_s < min_s)) {
            min_x = w3_x; 
            min_s = w3_s;   
            chord_side = 1;
            printf("je suis ici");
        }
        if (w4_x < min_x || (((w4_x - min_x)<1e-10) && w4_s <= min_s)) {
            printf("CHORD SIDE = %f\n", chord_side);
            min_x = w4_x; 
            min_s = w4_s;
            chord_side = -1;
        }
        printf("chord   side = %f\n", chord_side);
        if (chord_side == 1){
            BoundingBox_set_l(box, l);   // length of the box
            BoundingBox_set_w(box, w);   // width of the box
        } else {
            BoundingBox_set_l(box, w);   // length of the box
            BoundingBox_set_w(box, l);   // width of the box
        }
        BoundingBox_set_rC(box, sqrt((min_s-sR)*(min_s-sR) + (min_x-xR)*(min_x-xR))); // length of position vector to first vertex
    }
    printf("BoundingBox_get_rC(box) = %f\n", BoundingBox_get_rC(box));
    double rC = BoundingBox_get_rC(box);
    BoundingBox_set_sin_tC(box, (min_x-xR) / rC);   // angle of position vector to first vertex
    BoundingBox_set_cos_tC(box, (min_s-sR) / rC);
    printf("BoundingBox_get_sin_tC(box) = %f\n", BoundingBox_get_sin_tC(box));
    printf("BoundingBox_get_cos_tC(box) = %f\n", BoundingBox_get_cos_tC(box));
    BoundingBox_set_sin_tb(box, sin_t);           // orientation of the box (angle of length wrt horizontal)
    BoundingBox_set_cos_tb(box, cos_t);
    BoundingBox_set_proj_l(box, rC * (cos_t*s1/rC + sin_t*x1/rC)); // projection of position vector on length: rC * (cos_t*cos_tC + sin_t*sin_tC)
    BoundingBox_set_proj_w(box, rC * (cos_t*x1/rC - sin_t*s1/rC)); // projection of position vector on width:  rC * (cos_t*sin_tC - sin_t*cos_tC)
}


#endif /* XCOLL_GEOM_TRAJ_CIRCULAR_H */