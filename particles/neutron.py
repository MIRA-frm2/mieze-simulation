# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Neutron class."""

import numpy as np
import numpy.random as r
import os


cwd = os.getcwd()

c = 0.31225  # sqrt(1-polarisierung²)
x = c * r.rand()
z = np.sqrt(c ** 2 - x ** 2)

initial_polarisation = np.array([x, 0.95, z])


class Neutron:
    """Implements neutrons and its properties."""

    def __init__(self, speed, position, polarisation=initial_polarisation):
        """

        Parameters
        ----------
        speed: float
            The speed of the neutrons.
        position: np.array
            Array of the initial position of the neutron as (x, y, z).
        polarisation: np.array
            The array of the initial polarisation.
        """
        self.speed = speed
        self.polarisation = polarisation
        self.position = position

    def get_pol(self):
        """Return the polarisation of the neutron."""
        return self.polarisation

    def reset_pol(self):
        """Reset the polarisation to the initial value."""
        self.polarisation = initial_polarisation

    def set_position_x(self, x):
        """Set the neutron position."""
        self.position = (x, self.position[1], self.position[2])

    def get_position_x(self):
        """Return the neutron position."""
        return self.position[0]
