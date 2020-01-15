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
import random

from simulation.particles.neutron import Neutron

from utils.helper_functions import load_obj, find_nearest, rotate, get_phi
from simulation.experiments.mieze.parameters import angular_distribution_in_radians, speed_std


cwd = os.getcwd()


class NeutronBeam:
    """Implements neutrons and its properties."""

    def __init__(self, beamsize, incrementsize, speed, totalflightlength):
        self.neutrons = []
        self.number_of_neutrons = None
        self.speed = speed
        self.incrementsize = incrementsize
        self.totalflightlength = totalflightlength
        self.beamsize = beamsize

        self.polarisation = dict()

        self.b_map = None

        self.x_range = None
        self.y_range = None
        self.z_range = None

    def initialize_computational_space(self, **kwargs):
        """Initialize the 3d discretized computational space."""
        x_start = kwargs.pop('x_start', 0)
        x_end = kwargs.pop('x_end', 1)
        x_step = kwargs.pop('x_step', 0.1)

        y_start = kwargs.pop('y_start', 0)
        y_end = kwargs.pop('y_end', 1)
        z_start = kwargs.pop('z_start', 0)
        z_end = kwargs.pop('z_end', 1)
        yz_step = kwargs.pop('yz_step', 0.1)

        self.x_range = np.arange(x_start, x_end + x_step, x_step)
        self.y_range = np.arange(y_start, y_end + yz_step, yz_step)
        self.z_range = np.arange(z_start, z_end + yz_step, yz_step)

    def create_neutrons(self, number_of_neutrons, distribution=True):
        """Initialize the neutrons with a specific distribution.

        Parameters
        ----------
        number_of_neutrons: int
            Number of neutrons to be simulated.
        distribution: boolean, optional
            If True, then the neutrons are uniformly distributed.
            If False, then the neutrons are not spread and all start at 0.
            Defaults to True.
        """
        while len(self.neutrons) < number_of_neutrons:
            c = 0.31225  # sqrt(1-polarisierung²)
            x = c * r.rand()
            z = np.sqrt(c ** 2 - x ** 2)

            polarisation = np.array([x, 0.95, z])

            if distribution:

                pos_y = random.gauss(0, self.beamsize / 5)
                pos_z = random.gauss(0, self.beamsize / 5)

                # ToDo: cut to less digits for the position
                # pos_y = round(random.gauss(0, self.beamsize / 5), ndigits=-int(np.log10(self.incrementsize)))
                # pos_z = round(random.gauss(0, self.beamsize / 5), ndigits=-int(np.log10(self.incrementsize)))

                # ToDo: Randomized velocities
                # speed = random.gauss(self.speed, speed_std)
                #
                # radial_speed = speed * random.gauss(0, np.tan(angular_distribution_in_radians))
                # phi = get_phi(pos_y, pos_z)
                #
                # neutron_velocity = np.array([speed,
                #                              radial_speed * np.cos(phi),
                #                              radial_speed * np.sin(phi)])

                speed = random.gauss(self.speed, speed_std)

                radial_speed = speed * np.tan(angular_distribution_in_radians)
                phi = get_phi(pos_y, pos_z)

                neutron_velocity = np.array([speed,
                                             radial_speed * np.cos(phi),
                                             radial_speed * np.sin(phi)])

            else:
                pos_y = 0
                pos_z = 0

                neutron_velocity = np.array([self.speed, 0, 0])

            position = np.asarray([0, pos_y, pos_z])

            neutron = Neutron(polarisation=polarisation, position=position, velocity=neutron_velocity)

            self.neutrons.append(neutron)

    def _time_in_field(self, velocity):
        """Compute the time spent in the field."""
        return self.incrementsize / velocity

    @staticmethod
    def _omega(b):
        gamma = 1.83247172e4
        return b * gamma

    @staticmethod
    def _precession_angle(time, b):
        """Compute the precession angle."""
        gamma = 1.83247172e4
        return gamma * np.linalg.norm(b) * time

    def _polarisation_change(self, neutron, b, time):
        """Change the polarisation for the respective neutron.

        Parameters
        ----------
        neutron:
            Instance of Neutron class
        b: np.array
            Array of the magnetic field.

        Returns
        -------

        """
        phi = self._precession_angle(time, b)
        return rotate(vector=neutron.polarisation, phi=phi, axis=b)

    def load_magnetic_field(self):
        """Load the magnetic field data."""
        self.b_map = load_obj('../../data/data_magnetic_field')

    def compute_beam(self):
        """Compute the polarisation of the beam along the trajectory."""

        for j in self.x_range:
            for neutron in self.neutrons:

                # self.check_neutron_in_beam(neutron)

                time = self._time_in_field(velocity=neutron.speed)
                # print(j)
                neutron.set_position_x(j)
                neutron.compute_position_yz(time)
                # print(neutron.position)

                neutron.polarisation = self._polarisation_change(
                    neutron,
                    self.b_map[(find_nearest(self.x_range, neutron.position[0], index=False),
                                find_nearest(self.y_range, neutron.position[1], index=False),
                                find_nearest(self.z_range, neutron.position[2], index=False),
                                )],
                    time)

                neutron.trajectory.append(neutron.position)

            self.polarisation[tuple(self.get_neutron_position())] = self.get_pol()

            # print(self.get_pol())
            print(self.neutrons[0].trajectory)

    def check_neutron_in_beam(self, neutron):
        """Check if it is in the calculated beam profile (y,z plane)"""
        # if not (neutron.position[0], neutron.position[1], neutron.position[2]) in self.b_map:
        #     self.neutrons.remove(neutron)
        pass

    def get_pol(self):
        """Get the average polarisation for the beam."""
        xpol = np.average([neutron.polarisation[0] for neutron in self.neutrons])
        ypol = np.average([neutron.polarisation[1] for neutron in self.neutrons])
        zpol = np.average([neutron.polarisation[2] for neutron in self.neutrons])

        return xpol, ypol, zpol

    def reset_pol(self):
        """Reset the polarisation for each neutron to the initial polarisation."""
        for neutron in self.neutrons:
            neutron.reset_pol()

    def get_neutron_position(self, index=0):
        """Return the position of the neutron at index in the neutron list."""
        return self.neutrons[index].position
