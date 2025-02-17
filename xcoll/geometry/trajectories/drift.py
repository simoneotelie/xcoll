# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np

import xobjects as xo

from ...general import _pkg_root
from ..c_init import GeomCInit


class DriftTrajectory(xo.Struct):
    """A trajectory describing a drift.

    In parametrised form, it is given by:
        s(𝜆) = s0 + 𝜆 cos(𝜃0)     𝜆  = 0..∞
        x(𝜆) = x0 + 𝜆 sin(𝜃0)     𝜃0 = -π..π
        𝜃(𝜆) = 𝜃0

    where (s0, x0) and 𝜃0 are the initial particle coordinates resp. angle, and
    𝜆 represents the distance travelled along the direction 𝜃0.

    In practice, we do not provide 𝜃0 but xp = tan(𝜃0).
    """

    s0 = xo.Float64
    x0 = xo.Float64
    sin_t0 = xo.Float64
    cos_t0 = xo.Float64
    tan_t0 = xo.Float64

    _depends_on = [GeomCInit]
    _extra_c_sources = [_pkg_root / 'geometry' / 'trajectories' / 'drift.h']

    _kernels = {'set_params': xo.Kernel(
                                c_name='DriftTrajectory_set_params',
                                args=[xo.Arg(xo.ThisClass, name="traj"),
                                      xo.Arg(xo.Float64, name="s0"),
                                      xo.Arg(xo.Float64, name="x0"),
                                      xo.Arg(xo.Float64, name="xp")],
                                ret=None)}

    def __init__(self, *args, **kwargs):
        xp = kwargs.pop('xp', False)
        theta0 = kwargs.pop('theta0', False)
        super().__init__(*args, **kwargs)
        if xp is not False:
            self.xp = xp
        elif theta0 is not False:
            self.theta0 = theta0

    def __str__(self):
        return f"DriftTrajectory(s0={self.s0}, x0={self.x0}, xp={self.xp})"

    @property
    def xp(self):
        return self.tan_t0

    @property
    def theta0(self):
        return self.round(np.arctan2(self.sin_t0, self.cos_t0))
