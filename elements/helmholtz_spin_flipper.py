# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""An implementation of the Helmholtz Spin Flipper."""

from experiments.mieze.parameters import HelmholtzSpinFlipper_position_HSF1, unit

from elements.base import BasicElement
from elements.coils import RealCoil


class HelmholtzSpinFlipper(BasicElement):

    def __init__(self, position=HelmholtzSpinFlipper_position_HSF1, **kwargs):
        """Inherit init from base class and define additional parameters."""
        super(HelmholtzSpinFlipper, self).__init__()

        mid_pos = position

        length = kwargs.get('length',  0.01 * unit.m)  # length of each coil
        coil_distance = kwargs.get('coil_distance',  0.045 * unit.m)
        windings = kwargs.get('windings', 33)
        radius = kwargs.get('radius', 0.055 * unit.m)  # Positions.R_HSF # Radius of each coil
        current = kwargs.get('current', 1.6 * unit.m)

        pos1 = mid_pos - coil_distance / 2.0
        pos2 = mid_pos + coil_distance / 2.0

        self.coil1 = RealCoil(coil_mid_pos=pos1, length=length, windings=windings, r=radius, current=current)

        self.coil2 = RealCoil(coil_mid_pos=pos2, length=length, windings=windings, r=radius, current=current)

    def b_field(self, x, y, z):
        """Compute the magnetic field given the position."""
        b1 = self.coil2.b_field(x, y, z)
        b2 = self.coil2.b_field(x, y, z)
        return b1 + b2
