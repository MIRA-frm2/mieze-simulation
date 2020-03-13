# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019, 2020 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Neutron class."""

import numpy as np
import numpy.random as r
import os

from utils.helper_functions import rotate

cwd = os.getcwd()

c = 0.31225  # sqrt(1-polarisierung²)
x = c * r.rand()
z = np.sqrt(c ** 2 - x ** 2)

initial_polarisation = np.array([x, 0.95, z])


class Neutron:
    """Implements neutrons and its properties."""

    def __init__(self, velocity, position, polarisation=initial_polarisation):
        """

        Parameters
        ----------
        velocity: np.array
            The velocity vector of the neutron.
        position: np.array
            Array of the initial position of the neutron as (x, y, z).
        polarisation: np.array
            The array of the initial polarisation.
        """
        self.velocity = velocity
        self.polarisation = polarisation
        self.wavelength = 3956 / self.speed

        self.position = np.asarray(position)
        self.trajectory = list()

    @property
    def speed(self):
        """Return the speed from the velocity."""
        # return np.linalg.norm(self.velocity)
        return self.velocity[0]

    # Position related methods
    def set_position_x(self, x_val):
        """Set the neutron position."""
        self.position = np.array([x_val, self.position[1], self.position[2]])

    def update_position(self, time_increment):
        """Update the neutron position.

        Parameters
        ----------
        time_increment: float
            Time increment for the simulation.
        """
        self.position += np.array([time_increment * self.speed, 0, 0])

    def get_position_x(self):
        """Return the neutron position."""
        return self.position[0]

    def update_position_yz(self, time_increment):
        """Compute the position on the yz axes, based on the time spent in the cell on x axis.

        Parameters
        ----------
        time_increment: float
            Time spent in the cell on x axis.
        """
        self.position[1] += time_increment * self.velocity[1]
        self.position[2] += time_increment * self.velocity[2]

    # Polarisation related methods
    def get_pol(self):
        """Return the polarisation of the neutron."""
        return self.polarisation

    def reset_pol(self):
        """Reset the polarisation to the initial value."""
        self.polarisation = initial_polarisation

    @staticmethod
    def _precession_angle(time, b):
        """Compute the precession angle."""
        gamma = 1.83247172e4
        return gamma * np.linalg.norm(b) * time

    def compute_polarisation(self, magnetic_field_vector, time):
        """Change the neutron polarisation.

        Parameters
        ----------
        magnetic_field_vector: np.array
            Array of the magnetic field.
        time: float
            Time spent in this magnetic field cell.

        Returns
        -------

        """
        phi = self._precession_angle(time, b=magnetic_field_vector)
        self.polarisation = rotate(vector=self.polarisation, phi=phi, axis=magnetic_field_vector)
        return np.asarray(self.polarisation)
