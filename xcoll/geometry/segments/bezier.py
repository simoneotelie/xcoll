# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np

import xobjects as xo

from ...general import _pkg_root
from ..c_init import GeomCInit
from ..trajectories import DriftTrajectory


class BezierSegment(xo.Struct):
    """Bézier segment, defined by a start and end point P1 and P2, and two control points that define the curve"""
    s1 = xo.Float64
    x1 = xo.Float64
    s2 = xo.Float64
    x2 = xo.Float64
    cs1 = xo.Float64
    cx1 = xo.Float64
    cs2 = xo.Float64
    cx2 = xo.Float64

    _depends_on = [GeomCInit]
    _extra_c_sources = [_pkg_root / 'geometry' / 'segments' / 'bezier.h']

    max_crossings = {DriftTrajectory: 3}

    def __str__(self):
        return f"BezierSegment(({self.s1:.3}, {self.x1:.3})-c-({self.cs1:.3}, {self.cx1:.3}) -- " \
             + f"({self.cs2:.3}, {self.cx2:.3})-c-({self.s2:.3}, {self.x2:.3}))"

    def get_vertices(self):
        return (self.s1, self.x1), (self.s2, self.x2)

    def _translate_inplace(self, ds, dx):
        self.s1 += ds
        self.x1 += dx
        self.s2 += ds
        self.x2 += dx
        self.cs1 += ds
        self.cx1 += dx
        self.cs2 += ds
        self.cx2 += dx

    def _rotate_inplace(self, ps, px, angle):
        c = np.cos(angle)
        s = np.sin(angle)
        self._translate_inplace(-ps, -px)
        new_s1 = self.s1 * c - self.x1 * s
        new_x1 = self.s1 * s + self.x1 * c
        new_s2 = self.s2 * c - self.x2 * s
        new_x2 = self.s2 * s + self.x2 * c
        self.s1 = new_s1
        self.x1 = new_x1
        self.s2 = new_s2
        self.x2 = new_x2
        new_cs1 = self.cs1 * c - self.cx1 * s
        new_cx1 = self.cs1 * s + self.cx1 * c
        new_cs2 = self.cs2 * c - self.cx2 * s
        new_cx2 = self.cs2 * s + self.cx2 * c
        self.cs1 = new_cs1
        self.cx1 = new_cx1
        self.cs2 = new_cs2
        self.cx2 = new_cx2
        self._translate_inplace(ps, px)

