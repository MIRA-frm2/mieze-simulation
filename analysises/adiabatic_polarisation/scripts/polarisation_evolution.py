# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019, 2020 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Main script that computes the polarisation for a given setup."""

import imageio
import numpy as np

from simulation.beamline.beamline_mieze import MiezeBeamline
from analysises.neutron_polarisation_simulation.scripts.plotting_scripts import plot_polarisation_vector

from simulation.experiments.mieze.parameters import default_beam_grid, BEAM_PROPERTIES, \
    total_simulation_time, angular_distribution_in_radians, wavelength_min, wavelength_max
from utils.helper_functions import save_data_to_file


def main():
    """Main program that computes the neutron beam in the MIEZE setup."""

    # # Initialize an object from the MIEZE class
    # experiment = Mieze(spin_flipper_distance=ELEMENTS_POSITIONS_RELATIVE["spin_flipper_distance"],
    #                    coil_set_distance=ELEMENTS_POSITIONS_RELATIVE["coil_set_distance"],
    #                    save_individual_data_sets=False)
    #
    # # Create the components of the beamline with their parameters
    # experiment.create_setup()
    #
    # # Initialize the computational space (grid) and compute the magnetic field for it
    # experiment.initialize_computational_space(**default_beam_grid)
    # experiment.calculate_static_b_field()

    # Initialize the neutron beam
    simulation = MiezeBeamline(beamsize=BEAM_PROPERTIES['beamsize'],
                               speed=BEAM_PROPERTIES['speed'],
                               total_simulation_time=total_simulation_time)

    simulation.initialize_computational_space(**default_beam_grid)
    simulation.initialize_time_evolution_space()

    initial_polarisation = np.array([0, 1, 0])

    # Simulate the actual beam trajectory and the polarisation thereof
    images = list()
    time_factor = 1
    for t_j in np.linspace(0, time_factor * simulation.total_simulation_time,
                           num=time_factor * int(simulation.total_simulation_time / simulation.t_step)):
        # Show progress
        print(int(t_j/simulation.t_step))

        # Compute varying magnetic field
        simulation.load_magnetic_field(data_file_at_time_instance='../data/data_magnetic_field_mieze')
        # simulation.load_magnetic_field(data_file_at_time_instance='./../../../data/data_magnetic_field')

        # Create neutrons at each time step for the neutron beam
        simulation.create_neutrons(number_of_neutrons=1, distribution=False,
                                   polarisation=initial_polarisation)

        # Adjust the beam
        simulation.collimate_neutrons(max_angle=angular_distribution_in_radians)
        simulation.monochromate_neutrons(wavelength_min=wavelength_min, wavelength_max=wavelength_max)

        # Compute the neutrons movement and polarisation in the beam
        simulation.compute_beam()

        # Compute the average polarisation, at a fixed location (and time)
        simulation.compute_average_polarisation()
        if simulation.polarisation:
            print(f'pol:{list(simulation.polarisation.values())[-1]}')

        # Add image
        if simulation.polarisation:
            images.append(plot_polarisation_vector(simulation.polarisation, time=t_j, normalize=True, length=0.025))

    # Put all images together
    imageio.mimsave('../results/polarisation_vector.gif', images, fps=10)

    # save_data_to_file(simulation.polarisation, './data/data_polarisation.csv')


if __name__ == "__main__":
    main()
