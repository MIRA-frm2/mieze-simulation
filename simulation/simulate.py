# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019, 2020 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Main script that computes the magnetic field for a given experiment in a time dependent manner."""

import matplotlib.pyplot as plt
import numpy as np

from simulation.beamline.beam import NeutronBeam
from analysises.adiabatic_polarisation.scripts.adiabacity_parameter_plot import compute_polarisation
from analysises.neutron_polarisation_simulation.scripts.plotting_scripts import plot_polarisation_vector

from experiments.mieze.parameters import ELEMENTS_POSITIONS_RELATIVE
from simulation.beamline.beamline_properties import BEAM_PROPERTIES
from simulation.parameters_simulation import absolute_x_position, default_beam_grid, total_simulation_time


def simulate(experiment_class):
    """Main program that computes the neutron beam in the MIEZE experimental_setup."""

    # Initialize an object from the MIEZE class
    experiment = experiment_class(spin_flipper_distance=ELEMENTS_POSITIONS_RELATIVE["spin_flipper_distance"],
                                  coil_set_distance=ELEMENTS_POSITIONS_RELATIVE["coil_set_distance"],
                                  save_individual_data_sets=False)

    # Create the components of the beamline with their parameters
    experiment.create_setup()

    # Initialize the computational space (grid) and compute the magnetic field for it
    experiment.initialize_computational_space(**default_beam_grid)
    experiment.calculate_static_b_field()

    # Initialize the neutron beam
    simulation = NeutronBeam(beamsize=BEAM_PROPERTIES['beamsize'],
                             speed=BEAM_PROPERTIES['neutron_speed'],
                             total_simulation_time=total_simulation_time)

    simulation.initialize_computational_space(**default_beam_grid)
    simulation.initialize_time_evolution_space()

    polarisation_data = list()

    # Simulate the actual beam trajectory and the polarisation thereof

    for t_j in np.linspace(0, total_simulation_time, num=1):
        positions = list(absolute_x_position)
        for pos_x in absolute_x_position:

            # Show progress
            print(f'Computing for position: {pos_x}')

            # Compute varying magnetic field
            experiment.calculate_varying_magnetic_field(t_j)
            simulation.load_magnetic_field(b_map=experiment.b)
            # print(BEAM_PROPERTIES['initial_polarisation'])

            # Create neutrons at each time step for the neutron beam
            simulation.create_neutrons(number_of_neutrons=BEAM_PROPERTIES['number_of_neutrons'],
                                       distribution=False,
                                       polarisation=BEAM_PROPERTIES['initial_polarisation'])

            # Adjust the beam
            simulation.collimate_neutrons(max_angle=BEAM_PROPERTIES['angular_distribution'])
            simulation.monochromate_neutrons(wavelength_min=BEAM_PROPERTIES['wavelength_min'],
                                             wavelength_max=BEAM_PROPERTIES['wavelength_max'])

            # Compute the neutrons movement and polarisation in the beam
            simulation.compute_beam()

            # Compute the average polarisation, at a fixed location (and time)
            simulation.compute_average_polarisation()

            magnetic_field_vector = simulation.get_magnetic_field(np.asarray([pos_x, 0, 0]))
            if simulation.polarisation:
                polarisation_vector = list(simulation.polarisation.values())[-1]
                polarisation_value = compute_polarisation(polarisation_vector, magnetic_field_vector)

                polarisation_data.append(polarisation_value)
            else:
                positions.remove(pos_x)

        plot_polarisation_vector(polarisation_data=simulation.polarisation, save_image=False)

        plt.plot(positions, polarisation_data)
        plt.show()
