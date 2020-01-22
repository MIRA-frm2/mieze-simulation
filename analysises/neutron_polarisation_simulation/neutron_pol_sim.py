# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Main script that computes the flow of the particles through the magnetic field."""

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import axes3d, Axes3D  # <-- Note the capitalization!

from simulation.beamline.beamline_mieze import MiezeBeamline

from simulation.experiments.mieze.parameters import startpoint, beamend, beamsize, step_x, speed, \
    angular_distribution_in_radians, number_of_neutrons, wavelength_min, wavelength_max, default_beam_grid

from utils.helper_functions import save_data_to_file, read_data_from_file


class MyPlotter:

    def __init__(self, magnetic_field_data='../data/data.csv', polarisation_data='../data/polarisation_data.csv'):
        self.x_range, self.y_range, self.z_range, self.bx, self.by, self.bz = read_data_from_file(magnetic_field_data)
        self.x_range, self.y_range, self.z_range, self.px, self.py, self.pz = read_data_from_file(polarisation_data)

    def plot_field_3d(self, **kwargs):
        # 3d figure
        fig = plt.figure()
        ax = fig.gca(projection='3d')

        # Plot the magnetic field
        ax.quiver(self.x_range, self.y_range, self.z_range,
                  self.bx, self.by, self.bz,
                  color='b', **kwargs)
        ax.quiver(self.x_range, self.y_range, self.z_range,
                  self.px, self.py, self.pz,
                  color='g', **kwargs)

        plt.xlabel('x')
        plt.ylabel('y')

        plt.show()


def plot_polarisation_vector(polarisation_data='../../data/data_polarisation.csv', **kwargs):
    x_range, y_range, z_range, px, py, pz = read_data_from_file(polarisation_data)

    # 3d figure
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    # Plot the magnetic field
    ax.quiver(x_range[1::2], y_range[1::2], z_range[1::2],
              px[1::2], py[1::2], pz[1::2],
              color='g', **kwargs)

    plt.xlabel('x')
    plt.ylabel('y')

    plt.show()


def plot_neutron_trajectories(simulation):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for neutron in simulation.neutrons:
        for point in neutron.trajectory:
            ax.scatter(point[0], point[1], point[2], marker='o')


def main():
    simulation = MiezeBeamline(beamsize=beamsize,
                               incrementsize=step_x,
                               speed=speed,
                               totalflightlength=beamend)

    # Define computational space
    # grid_size = {'x_start': startpoint, 'x_end': beamend, 'x_step': step_x,
    #              'y_start': 0.0, 'y_end': 0.0, 'z_start': -0.0, 'z_end': 0.0,
    #              'yz_step':  (1.0 - -1.0) / 20}
    grid_size = default_beam_grid
    simulation.initialize_computational_space(**grid_size)

    # Load magnetic field values
    simulation.load_magnetic_field()

    # Initialize the neutrons and set their polarisation
    c = 0.31225  # sqrt(1-polarisierung²)
    x = c * np.random.rand()
    z = np.sqrt(c ** 2 - x ** 2)

    polarisation = np.array([x, 0.95, z])

    simulation.create_neutrons(number_of_neutrons=number_of_neutrons, distribution=False, polarisation=polarisation)
    simulation.reset_pol()

    # Adjust the beam
    simulation.collimate_neutrons(max_angle=angular_distribution_in_radians)
    simulation.monochromate_neutrons(wavelength_min=wavelength_min, wavelength_max=wavelength_max)

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

    plt.show()


if __name__ == "__main__":
    # main()
    plot_polarisation_vector(normalize=True, length=0.025)
