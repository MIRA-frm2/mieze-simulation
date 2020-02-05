# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019, 2020 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""An implementation of the polariser."""

import logging
import numpy as np

from simulation.elements.base import BasicElement
from experiments.mieze.parameters import POLARISATOR

from utils.physics_constants import MU_0, factor_T_to_G


# Create a custom logger
logger = logging.getLogger(__name__)


class Polariser(BasicElement):

    def __init__(self, position=(POLARISATOR, 0, 0), **kwargs):

        super(Polariser, self).__init__(position, name='Polariser')

        self.position = position

        self.c = kwargs.get('c', 0.05)  # shift polarisator position
        # magnetic dipol moment Polarisator;
        # both obtained by real measurements and data analysis
        # Pointing in the z upwards directions # Todo: check z or y
        self.m = kwargs.get('m', (0, 90, 0))

    def meta_data(self):
        return {"position": self.position, "c": self.c, "m": self.m}

    def _rectify_x_position(self, x):
        """Take into account the position and the shift on the x axis.

        The shift is that the magnetic field is measured/computed starting from the right side of the polarise.r
        """
        x -= self.position_x
        x += self.c
        return x

    def b_field(self, r: '(x, y, z)'):
        """Compute the magnetic field."""
        x, y, z = r
        if type(x) not in (int, float, np.float64):
            b_field = list()
            for item in x:
                b_field.append(np.abs(self.b_field([item, 0, 0])[1]))
            return np.asarray(b_field)
        else:
            x = self._rectify_x_position(x)

            r_vec = np.array((x, y, z))

            return self.b_field_theoretical(r_vec)

    def b_field_fitted(self, x_data, power, amplitude):
        """Compute the magnetic field as a power law, with parameters obtained from measured data."""
        x = self._rectify_x_position(np.asarray(x_data))
        return amplitude * np.power(x, -power) * 1e-4

    def b_field_theoretical(self, r_vec):
        """Compute the magnetic field from the equations of a magnetic dipole."""
        r = np.linalg.norm(r_vec)
        r_unit_vec = r_vec / r

        prefactor = MU_0 / (4 * np.pi)
        elem1 = 3 * r_unit_vec * np.asarray(self.m) * r_unit_vec
        elem2 = np.asarray(self.m)

        # Modelled as an ideal dipole magnetic field (note: it is divergent at 0)
        b = prefactor * (elem1 - elem2) / r**3
        return np.abs(b * factor_T_to_G)
