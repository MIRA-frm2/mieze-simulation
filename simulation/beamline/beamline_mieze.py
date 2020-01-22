# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2020 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Mieze beamline implementation."""

import numpy as np
import random


from .beam import NeutronBeam


from simulation.particles.neutron import Neutron
from simulation.experiments.mieze.parameters import angular_distribution_in_radians, speed_std

from utils.helper_functions import get_phi


class MiezeBeamline(NeutronBeam):

    def create_neutrons(self, distribution, number_of_neutrons, polarisation):

        """Initialize the neutrons with a specific distribution.

        Parameters
        ----------
        distribution: boolean, optional
            # ToDo: Requires more testing!
            If True, then the neutrons are uniformly distributed.
            If False, then the neutrons are not spread and all start at 0.
        number_of_neutrons: int
            Number of neutrons to be simulated.
        polarisation: np.array
            The initial polarisation of neutrons.
        """

        while len(self.neutrons) < number_of_neutrons:
            if distribution:

                # ToDo:
                # Make sure tha the magnetic field is computed at the computational grid required point.
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

                # Todo: Implement radial velocity such that some neutrons are still within the beamline at th end
                # radial_speed = speed * np.tan(angular_distribution_in_radians)
                radial_speed = 0
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
