# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Scripts used for plotting the results."""


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D  # <-- Note the capitalization!

from utils.helper_functions import read_data_from_file


def plot_polarisation_vector(polarisation_data='../../data/data_polarisation.csv', **kwargs):
    """Plot the polarisation vector along the beamline trajectory."""
    x_range, y_range, z_range, px, py, pz = read_data_from_file(polarisation_data)

    # 3d figure
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    # Plot the magnetic field
    ax.quiver(x_range, y_range, z_range,
              px, py, pz,
              color='g', **kwargs)

    plt.xlabel('x')
    plt.ylabel('y')

    plt.show()


def plot_neutron_trajectories(simulation):
    """Plot the trajectoies of the neutron in simulation."""
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for neutron in simulation.neutrons:
        for point in neutron.trajectory:
            ax.scatter(point[0], point[1], point[2], marker='o')
