# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""An implementation of the Helmholtz Spin Flipper."""

from experiments.mieze.parameters import HelmholtzSpinFlipper_position_HSF1
import numpy as np
from elements.coils import RealCoil


class HelmholtzSpinFlipper:
    MU_0 = 4e-7 * np.pi  # N/A, vacuum permeability
    R = 0.055  # Positions.R_HSF # Radius of each coil
    N = 33  # Windigs
    l = 0.01  # length of each coil
    d = 0.045  # coil distance

    def __init__(self):
        pass

    @classmethod
    def hsf(cls, x, y, z, mid_pos=HelmholtzSpinFlipper_position_HSF1, current=1.6):
        pos1 = mid_pos - cls.d/2.0
        pos2 = mid_pos + cls.d/2.0

        coil1 = RealCoil()
        coil1.create_element(coil_mid_pos=pos1, length=cls.l, windings=cls.N, r=cls.R)

        coil2 = RealCoil()
        coil2.create_element(coil_mid_pos=pos2, length=cls.l, windings=cls.N, r=cls.R)

        B1 = coil2.b_field(x, y, z)
        B2 = coil2.b_field(x, y, z)

        return B1 + B2
