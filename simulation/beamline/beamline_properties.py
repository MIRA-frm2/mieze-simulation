# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019, 2020 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""MIRA specific beamline properties."""

import numpy as np

# Beam properties
beamsize = 0.02
number_of_neutrons = 1


# Neutron wavelength and speed properties
wavelength = 4.3  # [Angstrom]

wavelength_min = 3.5
wavelength_max = 6

neutron_speed = 3956 / wavelength
speed_min = 3956 / wavelength_max
speed_max = 3956 / wavelength_min
speed_std = (speed_max - speed_min) / 4


c = 0.31225  # sqrt(1-polarisierung²)
x = c * np.random.rand()
z = np.sqrt(c ** 2 - x ** 2)

initial_polarisation = np.array([x, 0.95, z])

# Angular distribution
angular_distribution = 45  # minutes
angular_distribution_in_radians = angular_distribution * np.pi / (60 * 180)  # radians


BEAM_PROPERTIES = {
    'angular_distribution': angular_distribution_in_radians,
    'beamsize': beamsize,
    'initial_polarisation': initial_polarisation,
    'neutron_speed': neutron_speed,
    'number_of_neutrons': number_of_neutrons,
    'wavelength_min': wavelength_min,
    'wavelength_max': wavelength_max,
}

