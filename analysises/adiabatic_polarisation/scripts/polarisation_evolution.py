# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019, 2020 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Main script that computes the polarisation evolution for a given experimental_setup for testing."""

import imageio
import numpy as np

from simulation.beamline.beam import NeutronBeam
from simulation.beamline.beamline_properties import BEAM_PROPERTIES

from analysises.neutron_polarisation_simulation.scripts.plotting_scripts import plot_polarisation_vector
from simulation.parameters_simulation import default_beam_grid, total_simulation_time


def main():
    """Main program that computes the neutron beam in the MIEZE experimental_setup."""
    # Initialize the neutron beam
    simulation = NeutronBeam(beamsize=BEAM_PROPERTIES['beamsize'],
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
        simulation.collimate_neutrons(max_angle=BEAM_PROPERTIES['max_angle'])
        simulation.monochromate_neutrons(wavelength_min=BEAM_PROPERTIES['wavelength_min'],
                                         wavelength_max=BEAM_PROPERTIES['wavelength_max'])

        # Compute the neutrons movement and polarisation in the beam
        simulation.compute_beam()

        # Compute the average polarisation, at a fixed location (and time)
        simulation.compute_average_polarisation()

        if simulation.polarisation:
            print(f'pol:{list(simulation.polarisation.values())[-1]}')
            # Add image
            images.append(plot_polarisation_vector(simulation.polarisation, time=t_j))

    # Put all images together
    imageio.mimsave('../results/polarisation_vector.gif', images, fps=10)

    # save_data_to_file(simulation.polarisation, './data/data_polarisation.csv')


if __name__ == "__main__":
    main()
