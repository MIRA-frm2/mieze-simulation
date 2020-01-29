# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019, 2020 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Plotting scripts from data file."""

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D  # Needed for 3d plotting

import numpy as np

from utils.helper_functions import read_data_from_file, find_list_length_of_different_items


class Plotter:

    def __init__(self, filename='../../data/data_magnetic_field.csv'):
        self.x_range, self.y_range, self.z_range, self.bx, self.by, self.bz = read_data_from_file(filename)
        # self.preadjust_values()

    def preadjust_values(self):
        self.bx = 1e-1 * np.asarray(np.abs(self.bx))
        self.by = np.abs(self.by)

    def plot_field_1d_scalar(self, component=None, **kwargs):

        # fig = plt.figure()
        # ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])  # main axes
        if component == 'x':
            plt.plot(self.x_range, self.bx)
        elif component == 'y':
            plt.plot(self.x_range, self.by)
        elif component == 'z':
            plt.plot(self.x_range, self.bz)

        if kwargs:
            xlabel = kwargs.get('xlabel', None)
            ylabel = kwargs.get('ylabel', None)

            plt.xlabel(xlabel)
            plt.ylabel(ylabel)

        # if self.x_ticks and self.x_ticks_labels:
        #     ax.set_xticks(self.x_ticks)
        #     ax.set_xticklabels(self.x_ticks_labels)

        plt.show()

    def plot_field_2d_clr_map(self, plane='xy'):
        # fig = plt.figure()
        # ax = fig.gca()
        f = None
        # Plot the magnetic field
        if plane == 'xy':
            x = self.x_range
            y = self.y_range

            f = self.bz
            # fy = self.by
        elif plane == 'yz':
            x = self.y_range
            y = self.z_range

            f = self.bz
        elif plane == 'xz':
            x = self.x_range
            y = self.z_range

            # fx = self.bz
            # fy = self.bx
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

    def plot_field_3d(self, **kwargs):
        # 3d figure
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        # Plot the magnetic field
        ax.quiver(self.x_range, self.y_range, self.z_range,
                  self.bx, self.by, self.bz,
                  color='b', **kwargs)

        # plt.title('Magnetic field of a straight wire')
        plt.xlabel('x')
        plt.ylabel('y')

        plt.show()

    @staticmethod
    def extra_coil_check(b_extra, b_polariser):
        if b_extra.units != b_polariser.units:
            raise RuntimeError(f'B field units do not match:\n extra: {b_extra.units}\n'
                               f'polariser: {b_polariser.units}')

        if b_extra > b_polariser:
            pass
            # raise RuntimeError("Field of the extra coils too strong. Try resetting the extra coils.")


if __name__ == "__main__":
    plotter = Plotter()

    plotter.plot_field_1d_scalar(component='x')
    plotter.plot_field_1d_scalar(component='y')
    plotter.plot_field_1d_scalar(component='z')

    # plotter.plot_field_2d_clr_map(plane='yz')
    # plotter.plot_field_2d_vec_map(plane='xy')
    # plotter.plot_field_3d()

    # plotter.plot_adiabatic_check()
