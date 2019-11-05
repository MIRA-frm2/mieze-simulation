# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""An implementation of the Helmholtz Spin Flipper."""

from experiments.mieze.parameters import HelmholtzSpinFlipper_position_HSF1

from elements.base import BasicElement
from elements.coils import RealCoil

from utils.physics_constants import MU_0


class HelmholtzSpinFlipper(BasicElement):
    MU_0 = MU_0  # N/A, vacuum permeability

    def __init__(self):
        """Inherit init from base class and define additional parameters."""
        super(HelmholtzSpinFlipper, self).__init__()

        self.coil1 = None
        self.coil2 = None

    def create_element(self, position=HelmholtzSpinFlipper_position_HSF1, current=1.6, **kwargs):
        """Create the physical element."""
        mid_pos = position

        length = kwargs.get('length',  0.01)  # length of each coil
        coil_distance = kwargs.get('coil_distance',  0.045)
        windings = kwargs.get('windings', 33)
        radius = kwargs.get('radius', 0.055)  # Positions.R_HSF # Radius of each coil
        current = kwargs.get('current', 1.6)

        pos1 = mid_pos - coil_distance / 2.0
        pos2 = mid_pos + coil_distance / 2.0

        self.coil1 = RealCoil()
        self.coil1.create_element(coil_mid_pos=pos1, length=length, windings=windings, r=radius, current=current)

        self.coil2 = RealCoil()
        self.coil2.create_element(coil_mid_pos=pos2, length=length, windings=windings, r=radius, current=current)

    def hsf(self, x, y, z):
        """Compute the magnetic field given the position."""
        b1 = self.coil2.b_field(x, y, z)
        b2 = self.coil2.b_field(x, y, z)
        return b1 + b2
