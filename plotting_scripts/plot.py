# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Plotting scripts from data file."""

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D  # Needed for 3d plotting
import numpy as np
from scipy.interpolate import griddata

from utils.helper_functions import read_data_from_file, find_list_length_of_different_items


class Plotter:

    def __init__(self):
        self.x_range, self.y_range, self.z_range, self.bx, self.by, self.bz = read_data_from_file()

    def plot_field_1d_scalar(self, component=None):

        fig = plt.figure()
        # ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])  # main axes

        plt.plot(self.x_range, self.bx)

        # if self.x_ticks and self.x_ticks_labels:
        #     ax.set_xticks(self.x_ticks)
        #     ax.set_xticklabels(self.x_ticks_labels)

        plt.show()

    def plot_field_2d_clr_map(self, plane='xy'):
        fig = plt.figure()
        ax = fig.gca()
        # Plot the magnetic field
        if plane == 'xy':
            x = self.x_range
            y = self.y_range

            fx = self.bx
            fy = self.by
        elif plane == 'yz':
            x = self.y_range
            y = self.z_range

            f = self.bz
        elif plane == 'xz':
            x = self.x_range
            y = self.z_range

            fx = self.bz
            fy = self.bx
        else:
            x, y, fx, fy = None, None, None, None

        lenx = find_list_length_of_different_items(x)
        leny = find_list_length_of_different_items(y)

        # print(lenx, leny, len(fx))
        b = [leny * [0] for i in range(lenx)]

        for i in range(lenx):
            for j in range(leny):
                b[i][j] = f[j * lenx + i]

        plt.imshow(b, aspect='auto', extent=[min(x), max(x), min(y), max(y)])
        plt.colorbar()

        plt.show()

    def plot_field_2d_vec_map(self, plane='xy'):
        fig = plt.figure()
        ax = fig.gca()
        # Plot the magnetic field
        if plane == 'xy':
            x = self.x_range
            y = self.y_range

            fx = self.bx
            fy = self.by
        elif plane == 'yz':
            x = self.y_range
            y = self.z_range

            fx = self.by
            fy = self.bz
        elif plane == 'xz':
            x = self.x_range
            y = self.z_range

            fx = self.bz
            fy = self.bx
        else:
            x, y, fx, fy = None, None, None, None

        lenx = find_list_length_of_different_items(x)
        leny = find_list_length_of_different_items(y)

        # print(lenx, leny, len(fx))
        b = [leny * [0] for i in range(lenx)]

        for i in range(lenx):
            for j in range(leny):
                b[i][j] = fx[j * lenx + i] + fy[j * lenx + i]

        plt.imshow(b, aspect='auto', extent=[min(x), max(x), min(y), max(y)])
        plt.colorbar()

        for j in range(len(x)):
            f = np.sqrt(fx[j] ** 2 + fy[j] ** 2)
            ax.quiver(x[j], y[j], fx[j]/f, fy[j]/f, color='b', scale=100)

        # plt.title('Magnetic field of a straight wire')
        plt.xlabel('x')
        plt.ylabel('y')

        plt.show()

    def plot_field_3d(self):
        # 3d figure
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        # Plot the magnetic field
        ax.quiver(self.x_range, self.y_range, self.z_range,
                  self.bx, self.by, self.bz,
                  color='b')

        # plt.title('Magnetic field of a straight wire')
        plt.xlabel('x')
        plt.ylabel('y')

        plt.show()


if __name__ == "__main__":
    plotter = Plotter()

    # plotter.plot_field_1d_scalar()
    # plotter.plot_field_2d_clr_map(plane='yz')
    plotter.plot_field_2d_vec_map(plane='xy')
    # plotter.plot_field_3d()
