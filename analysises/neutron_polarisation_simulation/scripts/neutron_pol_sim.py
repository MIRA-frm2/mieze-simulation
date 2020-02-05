# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019, 2020 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Main script that computes the flow of the particles through the magnetic field."""

import numpy as np

from simulation.beamline.beam import NeutronBeam

from simulation.beamline.beamline_properties import BEAM_PROPERTIES
from simulation.parameters_simulation import beamend, default_beam_grid

from utils.helper_functions import save_data_to_file


def compute_neutron_beam():
    """Simulate the neutrons trajectories in the beamline."""
    total_simulation_time = beamend / BEAM_PROPERTIES['neutron_speed']
    simulation = NeutronBeam(beamsize=BEAM_PROPERTIES['beamsize'],
                             speed=BEAM_PROPERTIES['neutron_speed'],
                             total_simulation_time=total_simulation_time)

    # Define computational space
    # grid_size = {'x_start': startpoint, 'x_end': beamend, 'x_step': step_x,
    #              'y_start': 0.0, 'y_end': 0.0, 'z_start': -0.0, 'z_end': 0.0,
    #              'yz_step':  (1.0 - -1.0) / 20}
    grid_size = default_beam_grid
    simulation.initialize_computational_space(**grid_size)
    simulation.initialize_time_evolution_space()

    # Initialize the neutrons and set their polarisation
    c = 0.31225  # sqrt(1-polarisierung²)
    x = c * np.random.rand()
    z = np.sqrt(c ** 2 - x ** 2)

    polarisation = np.array([x, 0.95, z])

    simulation.create_neutrons(number_of_neutrons=BEAM_PROPERTIES['number_of_neutrons'],
                               distribution=False,
                               polarisation=polarisation)
    simulation.reset_pol()

    # Adjust the beam
    simulation.collimate_neutrons(max_angle=BEAM_PROPERTIES['angular_distribution'])
    simulation.monochromate_neutrons(wavelength_min=BEAM_PROPERTIES['wavelength_min'],
                                     wavelength_max=BEAM_PROPERTIES['wavelength_max'])

    # Simulate the actual beam trajectory and the polarisation thereof
    simulation.compute_beam()

    # Save data to file and plot
    save_data_to_file(simulation.polarisation, '../../data/data_polarisation')

    # plotter = MyPlotter(magnetic_field_data='../../data/data_magnetic_field.csv',
    #                     polarisation_data='../../data/data_polarisation.csv')
    # plotter.plot_field_3d(normalize=True, length=0.025)

    # plot_neutron_trajectories(simulation)
    # ax.set_ylim3d(-0.005, +0.005)
    # ax.set_zlim3d(-0.005, +0.005)

    # plt.show()


if __name__ == "__main__":
    compute_neutron_beam()
