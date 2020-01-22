# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019, 2020 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""General beamline class implementation."""

import numpy as np
import os
from abc import abstractmethod

from utils.helper_functions import load_obj, find_nearest, rotate


cwd = os.getcwd()


class NeutronBeam:
    """Implements neutrons and its properties."""

    def __init__(self, beamsize, speed):

        self.beamsize = beamsize
        self.speed = speed

        self.neutrons = []
        self.number_of_neutrons = None

        self.polarisation = dict()

        self.b_map = None

        self.x_range = None
        self.y_range = None
        self.z_range = None

        self.x_end = None

        self.x_step = None

        self.y_start = None
        self.y_end = None

        self.z_start = None
        self.z_end = None

        self.t_step = None
        self.t_end = None

    def initialize_computational_space(self, **kwargs):
        """Initialize the 3d discretized computational space.

        It should be identical (or at least similar) to the computed magnetic field values.
        """
        x_start = kwargs.pop('x_start', 0)
        self.x_end = kwargs.pop('x_end', 1)
        self.x_step = kwargs.pop('x_step', 0.1)

        self.y_start = kwargs.pop('y_start', 0)
        self.y_end = kwargs.pop('y_end', 1)
        self.z_start = kwargs.pop('z_start', 0)
        self.z_end = kwargs.pop('z_end', 1)

        yz_step = kwargs.pop('yz_step', 0.1)

        self.x_range = np.arange(x_start, self.x_end + self.x_step, self.x_step)
        self.y_range = np.arange(self.y_start, self.y_end + yz_step, yz_step)
        self.z_range = np.arange(self.z_start, self.z_end + yz_step, yz_step)

    def initialize_time_evolution_space(self):
        self.t_step = self.x_step / self.speed
        self.t_end = self.x_end / self.speed

        print(f'{self.t_step} {self.t_end}')

    @abstractmethod
    def create_neutrons(self, number_of_neutrons, polarisation, distribution):
        pass

    def _time_in_field(self, speed):
        """Compute the time spent in the field."""
        return self.x_step / speed

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

        for t_j in np.linspace(0, self.t_end, num=int(self.t_end / self.t_step)):
            for neutron in self.neutrons:

                self.check_neutron_in_beam(neutron)

                time = self._time_in_field(speed=neutron.speed)

                neutron.set_position_x(t_j * neutron.speed)
                neutron.compute_position_yz(time)

                magnetic_field_value = self.get_magnetic_field_value_at_neutron_position(neutron)

                neutron.polarisation = self._polarisation_change(
                    neutron,
                    magnetic_field_value,
                    time)

                neutron.trajectory.append(neutron.position)
            print(neutron.position)
            self.polarisation[tuple(neutron.position)] = self.get_pol()

    def get_magnetic_field_value_at_neutron_position(self, neutron):
        """Returns the magnetic field at the location of the magnetic field.

        Handles the case when the magnetic field does not contain a value at the required point.
        """
        try:
            return self.b_map[(find_nearest(self.x_range, neutron.position[0], index=False),
                               find_nearest(self.y_range, neutron.position[1], index=False),
                               find_nearest(self.z_range, neutron.position[2], index=False))]
        except KeyError:
            raise Exception(f'Could not find the magnetic field at the neutron position: {neutron.position}\n'
                            f'It is most probable that the magnetic field needs to be reevaluated.')

    def check_neutron_in_beam(self, neutron):
        """Check if neutron is in the calculated beam profile (y,z plane)"""
        if not (self.y_start <= neutron.position[1] <= self.y_end) \
                or not (self.z_start <= neutron.position[2] <= self.z_end):
            print("removed neutron")
            self.neutrons.remove(neutron)

    def get_pol(self):
        """Get the average polarisation for the beam."""
        # print([neutron.polarisation for neutron in self.neutrons])
        pol_x, pol_y, pol_z = 0, 0, 0
        n = len(self.neutrons)
        for neutron in self.neutrons:
            pol_x += neutron.polarisation[0]
            pol_y += neutron.polarisation[1]
            pol_z += neutron.polarisation[2]
        pol_x /= n
        pol_y /= n
        pol_z /= n

        # ToDo: compare which method is more computationally efficient
        # pol_x = np.average([neutron.polarisation[0] for neutron in self.neutrons])
        # pol_y = np.average([neutron.polarisation[1] for neutron in self.neutrons])
        # pol_z = np.average([neutron.polarisation[2] for neutron in self.neutrons])

        return pol_x, pol_y, pol_z

    def reset_pol(self):
        """Reset the polarisation for each neutron to the initial polarisation."""
        for neutron in self.neutrons:
            neutron.reset_pol()

    def get_neutron_position(self, index=0):
        """Return the position of the neutron at index in the neutron list."""
        return self.neutrons[index].position

    def collimate_neutrons(self, max_angle):
        """Apply a cut on the neutrons based on their angular distribution."""
        for neutron in self.neutrons:
            radial_speed = neutron.velocity[1] + neutron.velocity[2]
            angle = np.arctan(radial_speed/neutron.velocity[0])

            if angle > max_angle:
                self.neutrons.remove(neutron)

    def monochromate_neutrons(self, wavelength_min, wavelength_max):
        """Apply a cut on the neutrons based on their wavelength/speed distribution."""
        for neutron in self.neutrons:
            if neutron.wavelength > wavelength_max or neutron.wavelength < wavelength_min:
                self.neutrons.remove(neutron)
