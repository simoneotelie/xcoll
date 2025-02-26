// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_GEOM_SEG_LINE_H
#define XCOLL_GEOM_SEG_LINE_H


/*gpufun*/
int8_t LineSegment_bounded_below(LineSegment seg){
    return 1;  // closed segment
}

/*gpufun*/
int8_t LineSegment_bounded_above(LineSegment seg){
    return 1;  // closed segment
}

/*gpufun*/
double LineSegment_func_s(LineSegment seg, double t){
    double s1 = LineSegment_get_s1(seg);
    double s2 = LineSegment_get_s2(seg);
    return (1-t)*s1 + t*s2;
}

/*gpufun*/
double LineSegment_func_x(LineSegment seg, double t){
    double x1 = LineSegment_get_x1(seg);
    double x2 = LineSegment_get_x2(seg);
    return (1-t)*x1 + t*x2;
}

/*gpufun*/
double LineSegment_deriv_s(LineSegment seg, double t){
    UNUSED(t);
    return LineSegment_get_s2(seg) - LineSegment_get_s1(seg);
}

/*gpufun*/
double LineSegment_deriv_x(LineSegment seg, double t){
    UNUSED(t);
    return LineSegment_get_x2(seg) - LineSegment_get_x1(seg);
}

/*gpufun*/
void LineSegment_bounding_box_s(LineSegment seg, double t1, double t2, double extrema[2]){
    double s1 = LineSegment_func_s(seg, t1);
    double s2 = LineSegment_func_s(seg, t2);
    extrema[0] = MIN(s1, s2);
    extrema[1] = MAX(s1, s2);
}

/*gpufun*/
void LineSegment_bounding_box_x(LineSegment seg, double t1, double t2, double extrema[2]){
    double x1 = LineSegment_func_x(seg, t1);
    double x2 = LineSegment_func_x(seg, t2);
    extrema[0] = MIN(x1, x2);
    extrema[1] = MAX(x1, x2);
}





// /*gpufun*/
// void LineSegment_crossing_drift(LineSegment seg, int8_t* n_hit, double* s, double s0, double x0, double xm){
//     // Get segment data
//     double s1 = LineSegment_get_s1(seg);
//     double x1 = LineSegment_get_x1(seg);
//     double s2 = LineSegment_get_s2(seg);
//     double x2 = LineSegment_get_x2(seg);
//     double denom = (x2 - x1) - (s2 - s1)*xm;
//     if (fabs(denom) < XC_EPSILON){
//         // Trajectory is parallel to the segment
//         if (fabs((x0 - x1)/(s0 - s1) - xm) < XC_EPSILON){
//             // Trajectory overlaps with the segment
//             // TODO: This is situational; we should return s1 if get_s_first and current_s if after current_s
//             //       For now we hit twice (because we go nor IN nor OUT)
//             s[*n_hit] = s1;
//             (*n_hit)++;
//             s[*n_hit] = s2;
//             (*n_hit)++;
//         } else {
//             // No crossing
//             return;
//         }
//     } else {
//         double t = (x0 - x1 - (s0 - s1)*xm) / denom;
//         if (t >= 0 && t <= 1){
//             s[*n_hit] = s1*(1-t) + s2*t;
//             (*n_hit)++;
//         }
//     }
// }


// /*gpufun*/ 
// void LineSegment_crossing_mcs(LineSegment seg, int8_t* n_hit, double* s, const double* Ax, const double Xo, void* params){
//     return grid_search_and_newton();
// }
#endif /* XCOLL_GEOM_SEG_LINE_H */