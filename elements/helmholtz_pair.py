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
from elements.coils import Coil

from utils.physics_constants import factor_T_to_G


class HelmholtzPair(BasicElement):
    adjustment_factor = 1

    def __init__(self, position=(HelmholtzSpinFlipper_position_HSF1, 0, 0), **kwargs):
        """Inherit init from base class and define additional parameters."""
        super(HelmholtzPair, self).__init__(position, name='HelmHoltzSpinFlipper1')

        coil_type = kwargs.get('coil_type', Coil)

        width = kwargs.get('width',  0.01)  # width of each coil
        windings = kwargs.get('windings', 33)
        radius = kwargs.get('radius', 0.055)  # Positions.R_HSF # Radius of each coil
        current = kwargs.get('current', 1.6)

        pos1 = (self.position_x - radius / 2.0, 0, 0)
        pos2 = (self.position_x + radius / 2.0, 0, 0)

        self.coil1 = coil_type(name='HSF1',
                               position=pos1,
                               current=current,
                               length=width,
                               r_eff=radius,
                               windings=windings,
                               wire_d=0)

        self.coil2 = coil_type(name='HSF2',
                               position=pos2,
                               current=current,
                               length=width,
                               windings=windings,
                               r_eff=radius,
                               wire_d=0)

    def b_field(self, x, y, z):
        """Compute the magnetic field given the position."""
        b1 = self.coil2.b_field(x, y, z)
        b2 = self.coil2.b_field(x, y, z)
        return (b1 + b2) * self.adjustment_factor
