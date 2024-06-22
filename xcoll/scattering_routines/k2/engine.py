# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np
import os
from pathlib import Path
import shutil

import xtrack as xt


_FILE_NAME = "temp_k2_colldb.dat"

class K2Engine:
    # The engine is a singleton
    def __new__(cls, **kwargs):
        if not hasattr(cls, 'instance'):
            cls.instance = super().__new__(cls)
            cls.instance._initialised = False
        return cls.instance

    def __init__(self, **kwargs):
        if(self._initialised):
            return
        self._cwd = False
        self._old_cwd = False
        self._initialised = True
        self._warning_given = False
        self._file = None
        self._collimator_dict = {}
        self._capacity = np.int32(kwargs.get('_capacity', 50000))
        self.seed = kwargs.get('seed', None)

    def _warn(self):
        if not self.instance._warning_given:
            print("Warning: Failed to import pyK2 (did you compile?).\n"
                + "_K2collimators will be installed but are not trackable.")
            self._warning_given = True


    @classmethod
    def start(cls, line, *, seed=None, twiss=None, cwd=None, nemitt_x=None, nemitt_y=None):
        from ...beam_elements.k2 import _K2Collimator, _K2Crystal
        from .sixtrack_input import create_dat_file
        cls()
        this = cls.instance

        try:
            from .pyk2f import pyk2_init
        except ImportError:
           this._warn()
           return

        if line.particle_ref.mass0 != xt.Particles.reference_from_pdg_id('proton').mass0:
            raise ValueError("K2 only supports proton beams!")

        if seed is None:
            if this._seed is None:
                this._seed = np.int32(abs(np.random.randint(1, int(1.e16))))
        else:
            this._seed = np.int32(abs(seed))
            # Setting a seed here does not count as manually setting it
            this._seed_set_manually = False
        print(f"Using seed {this._seed}.")

        if cwd is not None:
            cwd = Path(cwd).expanduser().resolve()
            cwd.mkdir(parents=True, exist_ok=True)
            this._old_cwd = Path.cwd()
            os.chdir(cwd)
        else:
            cwd = Path.cwd()
        this._cwd = cwd

        elements, names = line.get_elements_of_type((_K2Collimator, _K2Crystal))

        if np.any([not hasattr(el, 'optics') or el.optics is None for el in elements]):
            raise ValueError("Not all collimators have optics assigned. Do this first")
        if nemitt_x is None:
            nemitt_x = np.unique([el.nemitt_x for el in elements])
            if len(nemitt_x) > 1:
                raise ValueError("Not all collimators have the same horizontal "
                               + "emittance. This is not supported.")
            nemitt_x = nemitt_x[0]
        if nemitt_y is None:
            nemitt_y = np.unique([el.nemitt_y for el in elements])
            if len(nemitt_y) > 1:
                raise ValueError("Not all collimators have the same vertical "
                               + "emittance. This is not supported.")
            nemitt_y = nemitt_y[0]

        if twiss is None:
            twiss = line.twiss()
        tw = twiss.rows[names]
        names = tw.name

        for i, name in enumerate(names):
            line[name]._k2_id = i + 1  # FORTRAN is 1-indexed

        this._file = create_dat_file(line=line, names=names, file=cwd / _FILE_NAME)
        assert this._file.exists()

        num_coll = np.int32(len(names))
        e0       = line.particle_ref.energy0[0]
        p0       = line.particle_ref.p0c[0]
        m0       = line.particle_ref.mass0
        beta0    = line.particle_ref.beta0[0]
        gamma0   = line.particle_ref.gamma0[0]

        pyk2_init(n_alloc=this._capacity, colldb_input_fname=this._file.name, random_generator_seed=this._seed, \
                  num_coll=num_coll, alphax=tw.alfx, alphay=tw.alfy, betax=tw.betx, betay=tw.bety, \
                  orbx=tw.x, orby=tw.y, orbxp=tw.px, orbyp=tw.py, nemitt_x=nemitt_x, nemitt_y=nemitt_y, \
                  e_ref=e0, p_ref=p0, m_ref=m0, beta_ref=beta0, gamma_ref=gamma0)

    @classmethod
    def stop(cls):
        cls()
        this = cls.instance
        this._cwd = None
        if this._old_cwd is not None:
            os.chdir(this._old_cwd)
            this._old_cwd = None
        if this._seed_set_manually:
            print("Warning: seed has not been reset. Do this manually before starting "
                  "the server again to avoid repeating the same random numbers (unless "
                  "this is intentional).")
        else:
            this._seed = None
        if this._file is not None:
            if this._file.exists():
                this._file.unlink()
            this._file = None

    @classmethod
    def is_running(cls):
        cls()
        this = cls.instance
        return this._file is not None or this._file.exists()


    @property
    def capacity(self):
        return self._capacity

    @capacity.setter
    def capacity(self, val):
        self._capacity = np.int32(val)


    @property
    def seed(self):
        return self._seed

    @seed.setter
    def seed(self, val):
        if val is None:
            self._seed_set_manually = False
        else:
            self._seed_set_manually = True
            val = np.int32(abs(val))
        self._seed = val

