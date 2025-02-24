# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np
import matplotlib.pyplot as plt

import xobjects as xo

from ..c_init import GeomCInit, PyMethod, XC_GEOM_EPSILON
from .drift import DriftTrajectory
from .circular import CircularTrajectory
from .mcs import MultipleCoulombTrajectory


all_trajectories = [DriftTrajectory, CircularTrajectory, MultipleCoulombTrajectory]


trajectory_methods = {
    'func_s': xo.Method(
        c_name=f"func_s",
        args=[xo.Arg(xo.Float64, name="l")],
        ret=xo.Arg(xo.Float64, name="s")),
    'func_x': xo.Method(
        c_name=f"func_x",
        args=[xo.Arg(xo.Float64, name="l")],
        ret=xo.Arg(xo.Float64, name="x")),
    'func_xp': xo.Method(
        c_name=f"func_xp",
        args=[xo.Arg(xo.Float64, name="l")],
        ret=xo.Arg(xo.Float64, name="theta")),
    'deriv_s': xo.Method(
        c_name=f"deriv_s",
        args=[xo.Arg(xo.Float64, name="l")],
        ret=xo.Arg(xo.Float64, name="s")),
    'deriv_x': xo.Method(
        c_name=f"deriv_x",
        args=[xo.Arg(xo.Float64, name="l")],
        ret=xo.Arg(xo.Float64, name="x")),
    'bounding_box_s': xo.Method(    # Gives the min/max s over the interval [t1, t2]
        c_name=f"bounding_box_s",
        args=[xo.Arg(xo.Float64, name="l1"),
              xo.Arg(xo.Float64, name="l2"),
              xo.Arg(xo.Float64, pointer=True, name="extrema"),
        ],
        ret=None),
    'bounding_box_x': xo.Method(    # Gives the min/max x over the interval [t1, t2]
        c_name=f"bounding_box_x",
        args=[xo.Arg(xo.Float64, name="l1"),
              xo.Arg(xo.Float64, name="l2"),
              xo.Arg(xo.Float64, pointer=True, name="extrema"),
        ],
        ret=None)
}


class LocalTrajectory(xo.UnionRef):
    """General trajectory, acting as a xobject-style parent class for all trajectory types"""
    _reftypes = all_trajectories
    _methods = list(trajectory_methods.values())

    def __init__(self, *args, **kwargs):
        raise ValueError("LocalTrajectory is an abstract class and should not be instantiated")

    @classmethod
    def from_dict(cls, dct, **kwargs):
        """Returns the correct trajectory object from a dictionary in the same style as a HybridClass"""
        this_dct = dct.copy()
        this_cls = this_dct.pop('__class__')
        class_found = False
        for cls in all_trajectories:
            if this_cls == cls.__name__:
                class_found = True
                break
        if not class_found:
            raise ValueError(f"Not a trajectory class: {this_cls}")
        return cls(**this_dct, **kwargs)


# Add kernels for func_ and deriv_ functions to all trajectories
def __getattr(self, attr):
    # Prepend the trajectory name to the kernel names to avoid duplication conflicts
    kernel_name = f"{self.__class__.__name__}_{attr}"
    if kernel_name in self._kernels:
        return PyMethod(kernel_name=kernel_name, element=self, element_name='traj')
    raise ValueError(f"Attribute {attr} not found in {self.__class__.__name__}")

for traj in all_trajectories:
    this_kernels = getattr(traj, '_kernels', {})
    _kernels = {key: xo.Kernel(c_name=f"{traj.__name__}_{val.c_name}",
                               ret=val.ret, args=[xo.Arg(xo.ThisClass, name="traj"), *val.args])
                for key, val in trajectory_methods.items()}
    this_kernels.update(_kernels)
    # Prepend the trajectory name to the kernel names to avoid duplication conflicts
    this_kernels = {f"{traj.__name__}_{key}": val for key, val in this_kernels.items()}
    traj._kernels = this_kernels
    traj.__getattr__ = __getattr
    traj._needs_compilation = True


# Define common methods for all trajectories
def plot(self, l1=0, l2=1):
    fig, ax  = plt.subplots()

    l_values = np.linspace(l1, l2, 100)
    s_values = np.array([self.func_s(l=l) for l in l_values])
    x_values = np.array([self.func_x(l=l) for l in l_values])

    # Plot the seg
    ax.plot(s_values, x_values, 'b-', label=f"{self.name} segment")

    # Plot the endpoints
    s_start, x_start = self.func_s(l=l1), self.func_x(l=l1)
    s_end, x_end = self.func_s(l=l2), self.func_x(l=l2)
    ax.plot([s_start, s_end], [x_start, x_end], 'go', label='Endpoints')

    # Get and plot the bounding box
    extrema_s = np.zeros(2)
    extrema_x = np.zeros(2)
    self.bounding_box_s(l1=l1, l2=l2, extrema=extrema_s)
    self.bounding_box_x(l1=l1, l2=l2, extrema=extrema_x)
    ax.plot([extrema_s[0], extrema_s[1], extrema_s[1], extrema_s[0], extrema_s[0]],
            [extrema_x[0], extrema_x[0], extrema_x[1], extrema_x[1], extrema_x[0]], 'k--', label='Bounding Box')

    # Set plot labels and show the plot
    ax.set_xlabel('s')
    ax.set_ylabel('x')
    ax.legend()

    return ax.figure, ax

def __eq(self, other):
    """Check if two objects are equal"""
    return self.to_dict() == other.to_dict()

def __repr(self):
    """Return repr(self)."""
    return f"<{str(self)} at {hex(id(self))}>"

def to_dict(self):
    """Returns a dictionary in the same style as a HybridClass"""
    return {'__class__': self.__class__.__name__, **self._to_json()}

@classmethod
def from_dict(cls, dct, **kwargs):
    """Returns the object from a dictionary in the same style as a HybridClass"""
    this_dct = dct.copy()
    this_cls = this_dct.pop('__class__')
    if this_cls != cls.__name__:
        raise ValueError(f"Expected class {cls.__name__}, got {this_cls}")
    return cls(**this_dct, **kwargs)

def __copy(self):
    """Returns a copy of the object"""
    return self.from_dict(self.to_dict())

def __round(self, val):
    """Built-in to provide rounding to Xcoll precision"""
    return round(val, -int(np.log10(XC_GEOM_EPSILON)))

for traj in all_trajectories:
    traj.name = traj.__name__.lower()[:-10]
    traj.__eq__ = __eq
    if not '__repr__' in traj.__dict__:
        traj.__repr__ = __repr
    if not '__str__' in traj.__dict__:
        traj.__str__ = traj__repr__
    traj.to_dict = to_dict
    traj.from_dict = from_dict
    traj.copy = __copy
    traj.round = __round
    traj.plot = plot


# Sanity check to assert all trajectories have C code for func_ and deriv_ functions
def assert_trajectory_sources(tra):
    assert traj in all_trajectories
    def _check_source(header, traj):
        header_found = False
        for src in traj._extra_c_sources:
            if isinstance(src, str):
                if header in src:
                    header_found = True
                    break
            else:
                with open(src) as f:
                    if header in f.read():
                        header_found = True
                        break
        if not header_found:
            raise SystemError(f"Missing or corrupt C function in trajectory {traj.__name__}:  {header}.")
    name = traj.__name__
    for func in ['func_s', 'func_x', 'func_xp', 'deriv_s', 'deriv_x']:
        header = f"/*gpufun*/\ndouble {name}_{func}({name} traj, double l)"
        _check_source(header, traj)
    for func in ['bounding_box_s', 'bounding_box_x']:
        header = f"/*gpufun*/\nvoid {name}_{func}({name} traj, double l1, double l2, double extrema[2])"
        _check_source(header, traj)
for traj in all_trajectories:
    assert_trajectory_sources(traj)




# OLD ========================================================================================================
from ..c_init import xo_to_ctypes

args_cross_h = [
    xo.Arg(xo.Int8,    pointer=True, name="n_hit"),
    xo.Arg(xo.Float64, pointer=True, name="s"),
]
args_cross_v = [
    xo.Arg(xo.Float64, pointer=True, name="restrict_s"),
]
args_vlimit = [
    xo.Arg(xo.Float64, pointer=False, name="ymin"),
    xo.Arg(xo.Float64, pointer=False, name="ymax")
]
args_length = [
    xo.Arg(xo.Float64, pointer=False, name="s1"),
    xo.Arg(xo.Float64, pointer=False, name="s2")
]




# # Sanity check to assert all trajectory types have a vlimit function and a length function
# def assert_trajectory_sources(tra):
#     assert tra in all_trajectories
#     header = f"/*gpufun*/\ndouble {tra.__name__}_func({xo_to_ctypes(args_cross_v)}, {xo_to_ctypes(tra.args_hv)}, " \
#             + f"{xo_to_ctypes(tra.args_v)}, {xo_to_ctypes(args_vlimit)})"
#     header_found = False
#     for src in tra._extra_c_sources:
#         if isinstance(src, str):
#             if header in src:
#                 header_found = True
#                 break
#         else:
#             with open(src) as f:
#                 if header in f.read():
#                     header_found = True
#                     break
#     if not header_found:
#         raise SystemError(f"Missing or corrupt C function for {tra.__name__}.")

#     header = f"/*gpufun*/\nint8_t {tra.__name__}_vlimit({xo_to_ctypes(args_cross_v)}, {xo_to_ctypes(tra.args_hv)}, " \
#             + f"{xo_to_ctypes(tra.args_v)}, {xo_to_ctypes(args_vlimit)})"
#     header_found = False
#     for src in tra._extra_c_sources:
#         if isinstance(src, str):
#             if header in src:
#                 header_found = True
#                 break
#         else:
#             with open(src) as f:
#                 if header in f.read():
#                     header_found = True
#                     break
#     if not header_found:
#         raise SystemError(f"Missing or corrupt C vlimit function for {tra.__name__}.")

#     header = f"/*gpufun*/\ndouble {tra.__name__}_length({xo_to_ctypes(tra.args_hv)}, {xo_to_ctypes(tra.args_h)}, " \
#                 f"{xo_to_ctypes(tra.args_v)}, {xo_to_ctypes(args_length)})"
#     header_found = False
#     for src in tra._extra_c_sources:
#         if isinstance(src, str):
#             if header in src:
#                 header_found = True
#                 break
#         else:
#             with open(src) as f:
#                 if header in f.read():
#                     header_found = True
#                     break
#     if not header_found:
#         raise SystemError(f"Missing or corrupt C length function for {tra.__name__}.")
