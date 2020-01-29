# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019, 2020 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""An implementation of the Spin Flipper."""

from simulation.experiments.mieze.parameters import SpinFlipper_position1
from numpy import arctan

from simulation.elements.coils import RectangularCoil

from utils.physics_constants import MU_0, pi


class SpinFlipper(RectangularCoil):
    """Class that implements a SpinFlipper coil."""

    class_name = 'SpinFlipper'

    def __init__(self, name='SpinFlipper', position=(SpinFlipper_position1, 0, 0), **kwargs):
        super(SpinFlipper, self).__init__(position, name=name, **kwargs)

        self.thickness = kwargs.get('thickness', 1e-2)  # [m], thickness of the coil

    def meta_data(self):
        """Return meta data."""
        return {"position": self.position_x, "coil_type": "RectangularCoil"}

    def _sf_th(self, current, position):
        """Compute spin flipper approximation."""
        x0 = self.length / 2.0
        n = self.windings / self.length

        if position == 0:
            angle = pi / 2
        else:
            angle = arctan(x0 / position)

        b_eff = n * MU_0 * current / pi * angle
        return b_eff
    
    def b_field_approx(self, position):
        """Compute an approximation of the magnetic field along the beam axis.

        Parameters
        ----------
        position: ndarray, list
            Position along the beam axis.

        Returns
        -------
        out: ndarray, list
            Magnetic field vector
        """
        x = position[0]
        y = x - self.position_x

        b1 = self._sf_th(self.current, y + self.thickness / 2.0)
        b2 = -self._sf_th(self.current, y - self.thickness / 2.0)
        return [0, (b1 + b2) * 1e4, 0]
