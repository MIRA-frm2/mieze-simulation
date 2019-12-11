# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Main script that computes the flow of the particles through the magnetic field."""

# import numpy as np

from particles.beam import NeutronBeam

from experiments.mieze.parameters import I_hsf1, step_x, startpoint, beamend, beamsize, step_x, velocity
from plotting_scripts.plot import Plotter

from utils.helper_functions import save_data_to_file
from utils.physics_constants import earth_field


# def b_function(vec):
#     x = vec[0]
#     y = vec[1]
#     z = vec[2]
#
#     flipper = HelmholtzSpinFlipper()
#     flipper.create_element(position=0.1, current=I_hsf1)
#
#     return flipper.hsf(x, y, z) + Polariser.B_field(x, y, z) + earth_field
#

def main():
    simulation = NeutronBeam(beamsize=beamsize,
                             incrementsize=step_x,
                             number_of_neutrons=1,
                             velocity=velocity,
                             totalflightlength=beamend)

    grid_size = {'x_start': startpoint, 'x_end': beamend, 'x_step': step_x,
                 'y_start': 0.0, 'y_end': 0.0, 'z_start': -0.0, 'z_end': 0.0,
                 'yz_step':  (1.0 - -1.0) / 20}

    simulation.initialize_computational_space(**grid_size)

    simulation.create_neutrons()
    simulation.reset_pol()

    print(simulation.get_pol())
    print(simulation.get_neutron_position())

    simulation.compute_beam()

    save_data_to_file(simulation.polarisation, '../data/data_polarisation')

    plotter = Plotter('../data/data_polarisation.csv')
    plotter.plot_field_3d(normalize=True, length=0.05)


if __name__ == "__main__":
    main()
