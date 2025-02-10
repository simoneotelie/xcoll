// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_GEOM_SEG_HALFOPENLINE_H
#define XCOLL_GEOM_SEG_HALFOPENLINE_H
#define XC_HALFOPENLINE_CROSSINGS 2



/*gpufun*/
double HalfOpenLineSegment_func(HalfOpenLineSegment seg, double s){
    // MCS trajectory form PDG rewritted in terms of A, B and s/Xo. 
    const double Ax = HalfOpenLine_get_Ax(seg);
    const double Xo = HalfOpenLine_get_Xo(seg);
    double s2       = HalfOpenLine_get_s2(seg);
    double x2       = HalfOpenLine_get_x2(seg);
    double s1       = HalfOpenLine_get_s1(seg);
    double x1       = HalfOpenLine_get_x1(seg);
    return x2 + (x2 - x1) / (s2 - s1) * (s - s1);
}

/*gpufun*/
double HalfOpenLineSegment_deriv(HalfOpenLineSegment seg, double s){
    // MCS trajectory derivative wrt s
    const double Ax = HalfOpenLine_get_Ax(seg);
    const double Xo = HalfOpenLine_get_Xo(seg);
    double s2       = HalfOpenLine_get_s2(seg);
    double x2       = HalfOpenLine_get_x2(seg);
    double s1       = HalfOpenLine_get_s1(seg);
    double x1       = HalfOpenLine_get_x1(seg);
    return (x2 - x1) / (s2 - s1);
}


/*gpufun*/
void HalfOpenLineSegment_crossing_drift(HalfOpenLineSegment seg, int8_t* n_hit, double* s, double s0, double x0, double xm){
    // Get segment data
    double s1 = HalfOpenLineSegment_get_s(seg);
    double x1 = HalfOpenLineSegment_get_x(seg);
    double s2 = s1 + cos(HalfOpenLineSegment_get_t(seg));
    double x2 = x1 + sin(HalfOpenLineSegment_get_t(seg));
    double denom = (x2 - x1) - (s2 - s1)*xm;
    if (fabs(denom) < XC_EPSILON){
        // Trajectory is parallel to the segment
        if (fabs((x0 - x1)/(s0 - s1) - xm) < XC_EPSILON){
            // Trajectory overlaps with the segment
            // TODO: This is situational; we should return s1 if get_s_first and current_s if after current_s
            //       For now we hit twice (because we go nor IN nor OUT)
            s[*n_hit] = s1;
            (*n_hit)++;
            s[*n_hit] = s1;
            (*n_hit)++;
        } else {
            // No hit
            return;
        }
    } else {
        double t = (x0 - x1 - (s0 - s1)*xm) / denom;
        if (t >= 0){  // We do not check for t<=1 as it is a half-open segment
            s[*n_hit] = s1*(1-t) + s2*t;
            (*n_hit)++;
        }
    }
}

/*gpufun*/
void HalfOpenLineSegment_crossing_mcs(HalfOpenLineSegment seg, int8_t* n_hit, double* s, const double* Ax, const double Xo){
   return grid_search_and_newton()
}  
#endif /* XCOLL_GEOM_SEG_HALFOPENLINE_H */