# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np

import xobjects as xo
from ..c_init import GeomCInit

# Bézier segment, defined by a start and end point P1 and P2, and two control points that define the curve
class BezierSegment(xo.Struct):
    s1 = xo.Float64
    x1 = xo.Float64
    s2 = xo.Float64
    x2 = xo.Float64
    cs1 = xo.Float64
    cx1 = xo.Float64
    cs2 = xo.Float64
    cx2 = xo.Float64

    _depends_on = [GeomCInit]
    _extra_c_sources = [
        """
/*gpufun*/
void _hit_s_bezier(BezierSegment seg, double t, double multiplicity, int8_t* n_hit, double* s){
    double s1  = BezierSegment_get_s1(seg);
    double s2  = BezierSegment_get_s2(seg);
    double cs1 = BezierSegment_get_cs1(seg);
    double cs2 = BezierSegment_get_cs2(seg);
    double new_s = (1-t)*(1-t)*(1-t)*s1 + 3*(1-t)*(1-t)*t*cs1 + 3*(1-t)*t*t*cs2 + t*t*t*s2;
    for (int8_t i = 0; i < multiplicity; i++) {
        s[*n_hit] = new_s;
        (*n_hit)++;
    }
}

/*gpufun*/
void BezierSegment_crossing_drift(BezierSegment seg, int8_t* n_hit, double* s, double s0, double x0, double m){
    // Get segment data
    double s1  = BezierSegment_get_s1(seg);
    double x1  = BezierSegment_get_x1(seg);
    double s2  = BezierSegment_get_s2(seg);
    double x2  = BezierSegment_get_x2(seg);
    double cs1 = BezierSegment_get_cs1(seg);
    double cx1 = BezierSegment_get_cx1(seg);
    double cs2 = BezierSegment_get_cs2(seg);
    double cx2 = BezierSegment_get_cx2(seg);
    // The Bézier curve is defined by the parametric equations (with t in [0, 1]):
    // s(t) = (1-t)^3*s1 + 3(1-t)^2*t*cs1 + 3(1-t)*t^2*cs2 + t^3*s2
    // x(t) = (1-t)^3*x1 + 3(1-t)^2*t*cx1 + 3(1-t)*t^2*cx2 + t^3*x2
    // Plug the parametric eqs into the drift trajectory x(t) = m*(s(t) - s0) + x0 and solve for t
    // The solutions for t (which we get by Cardano's method) are valid if in [0, 1]
    double a = (m*s1 - x1) - (m*s2 - x2) - 3*(m*cs1 - cx1) + 3*(m*cs2 - cx2);
    double b = 6*(m*cs1 - cx1) - 3*(m*cs2 - cx2) - 3*(m*s1 - x1);
    double c = 3*(m*s1 - x1) - 3*(m*cs1 - cx1);
    double d = (m*s0 - x0) - (m*s1 - x1);
    double t;
    // Edge cases
    if (fabs(a) < XC_EPSILON){
        if (fabs(b) < XC_EPSILON){
            if (fabs(c) < XC_EPSILON){
                if (fabs(d) < XC_EPSILON){
                    // The trajectory is on the Bézier curve
                    // TODO: This cannot happen because we don't have a cubic trajectory.
                    //       Adapt if these ever would be implemented.
                    return;
                } else {
                    // No solutions
                    return;
                }
            } else {
                // This is a linear equation
                t = -d/c;
                if (0 <= t && t <= 1){
                    _hit_s_bezier(seg, t, 1, n_hit, s);
                }
            }
        } else {
            // This is a quadratic equation
            double disc = c*c - 4*b*d;
            if (disc < 0){
                // No solutions
                return;
            }
            for (int8_t i = 0; i < 2; i++) {
                double sgnD = i*2-1; // negative and positive solutions; if multiplicity 2, we add the same solution twice
                t = (-c + sgnD*sqrt(fabs(disc)))/(2*b);
                if (0 <= t && t <= 1){
                    _hit_s_bezier(seg, t, 1, n_hit, s);
                }
            }
        }
    } else {
        // Full cubic equation. Coefficients for the depressed cubic t^3 + p*t + q = 0:
        double p = (3*a*c - b*b)/(3*a*a);
        double q = (2*b*b*b - 9*a*b*c + 27*a*a*d)/(27*a*a*a);
        double disc = -p*p*p/27 - q*q/4;  // This is the discriminant of the depressed cubic but divided by (4*27)
        if (fabs(disc) < XC_EPSILON){
            if (fabs(p) < XC_EPSILON){
                // One real root with multiplicity 3
                t = -b/(3*a);
                if (0 <= t && t <= 1){
                    _hit_s_bezier(seg, t, 3, n_hit, s);
                }
            } else {
                // Two real roots (one simple and one with multiplicity 2)
                t = 3*q/p - b/(3*a);
                if (0 <= t && t <= 1){
                    _hit_s_bezier(seg, t, 1, n_hit, s);
                }
                t = -3*q/(2*p) - b/(3*a);
                if (0 <= t && t <= 1){
                    _hit_s_bezier(seg, t, 2, n_hit, s);
                }
            }
        } else if (disc < 0){
            // One real root
            t = cbrt(-q/2 + sqrt(fabs(disc))) + cbrt(-q/2 - sqrt(fabs(disc))) - b/(3*a);
            if (0 <= t && t <= 1){
                _hit_s_bezier(seg, t, 1, n_hit, s);
            }
        } else {
            // Three real roots
            double phi = acos(3*q/(2*p)*sqrt(fabs(3/p)));
            t = 2*sqrt(fabs(p/3))*cos(phi/3) - b/(3*a);
            if (0 <= t && t <= 1){
                _hit_s_bezier(seg, t, 1, n_hit, s);
            }
            t = 2*sqrt(fabs(p/3))*cos((phi + 2*M_PI)/3) - b/(3*a);
            if (0 <= t && t <= 1){
                _hit_s_bezier(seg, t, 1, n_hit, s);
            }
            t = 2*sqrt(fabs(p/3))*cos((phi + 4*M_PI)/3) - b/(3*a);
            if (0 <= t && t <= 1){
                _hit_s_bezier(seg, t, 1, n_hit, s);
            }
        }
    }
}
"""
    ]

    def evaluate(self, t):
        s1  = self.s1
        x1  = self.x1
        s2  = self.s2
        x2  = self.x2
        cs1 = self.cs1
        cx1 = self.cx1
        cs2 = self.cs2
        cx2 = self.cx2
        t = np.array(t)
        mask = (t >= 0) & (t <= 1)
        return (1-t[mask])**3 * s1 + 3*t[mask]*(1-t[mask])**2 * cs1 + 3*(1-t[mask])*t[mask]**2 * cs2 + t[mask]**3 * s2, \
               (1-t[mask])**3 * x1 + 3*t[mask]*(1-t[mask])**2 * cx1 + 3*(1-t[mask])*t[mask]**2 * cx2 + t[mask]**3 * x2