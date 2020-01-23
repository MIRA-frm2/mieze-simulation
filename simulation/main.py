# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Main script that computes the magnetic field for a given experiment."""

import imageio
import numpy as np
import matplotlib.pyplot as plt

from simulation.experiments.mieze.main_mieze import Mieze
from simulation.beamline.beamline_mieze import MiezeBeamline
from analysises.neutron_polarisation_simulation.scripts.plotting_scripts import plot_polarisation_vector

from simulation.experiments.mieze.parameters import ELEMENTS_POSITIONS_RELATIVE, default_beam_grid, beamsize, speed, \
    total_simulation_time, initial_polarisation, number_of_neutrons, angular_distribution_in_radians, wavelength_min, \
    wavelength_max


def main():
    # Initialize an object from the MIEZE class
    experiment = Mieze(spin_flipper_distance=ELEMENTS_POSITIONS_RELATIVE["spin_flipper_distance"],
                       coil_set_distance=ELEMENTS_POSITIONS_RELATIVE["coil_set_distance"],
                       save_individual_data_sets=False)

    # Create the components of the beamline with their parameters
    experiment.create_setup()

    # Initialize the computational space (grid) and compute the magnetic field for it
    experiment.initialize_computational_space(**default_beam_grid)
    experiment.calculate_b_field()

    # Initialize the neutron beam
    simulation = MiezeBeamline(beamsize=beamsize, speed=speed, total_simulation_time=total_simulation_time)

    simulation.initialize_computational_space(**default_beam_grid)
    simulation.initialize_time_evolution_space()

    # Initialize the neutrons and set their polarisation
    simulation.create_neutrons(number_of_neutrons=number_of_neutrons, distribution=False,
                               polarisation=initial_polarisation)
    simulation.reset_pol()

    # Adjust the beam
    simulation.collimate_neutrons(max_angle=angular_distribution_in_radians)
    simulation.monochromate_neutrons(wavelength_min=wavelength_min, wavelength_max=wavelength_max)

    # Simulate the actual beam trajectory and the polarisation thereof
    images = list()
    i = 0
    for t_j in np.linspace(0, simulation.total_simulation_time, num=int(simulation.total_simulation_time / simulation.t_step)):
        simulation.load_magnetic_field(t_j, b_map=experiment.b)
        simulation.compute_beam(t_j)

        print(i)
        images.append(plot_polarisation_vector(simulation.polarisation, normalize=True, length=0.025))
        i += 1

    imageio.mimsave('./test.gif', images, fps=10)


if __name__ == "__main__":
    main()
