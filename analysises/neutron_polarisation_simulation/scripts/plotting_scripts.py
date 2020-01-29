# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019, 2020 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Scripts used for plotting the results."""


import matplotlib.pyplot as plt
# Needed implicitly for the 3d plots, although not explicitly used
from mpl_toolkits.mplot3d import axes3d, Axes3D
import numpy as np

from utils.helper_functions import read_data_from_file


def plot_polarisation_vector(polarisation_data=None, time=None,
                             polarisation_data_file='../../data/data_polarisation.csv', **kwargs):
    """Plot the polarisation vector along the beamline trajectory."""

    # 3d figure
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    plt.title(f'Time: {time}')
    plt.xlabel('x [m]')
    plt.ylabel('y')

    # Plot the magnetic field
    if polarisation_data:
        for keys, values in polarisation_data.items():
            ax.quiver(keys[0], keys[1], keys[2],
                      values[0], values[1], values[2],
                      color='g', **kwargs)

        # Used to return the plot as an image rray
        fig.canvas.draw()       # draw the canvas, cache the renderer
        image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
        image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))

        plt.close(fig)

        return image

    elif polarisation_data_file:
        x_range, y_range, z_range, px, py, pz = read_data_from_file(polarisation_data)
        ax.quiver(x_range, y_range, z_range,
                  px, py, pz,
                  color='g', **kwargs)
        plt.show()


def plot_neutron_trajectories(simulation):
    """Plot the trajectoies of the neutron in simulation."""
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for neutron in simulation.neutrons:
        for point in neutron.trajectory:
            ax.scatter(point[0], point[1], point[2], marker='o')
