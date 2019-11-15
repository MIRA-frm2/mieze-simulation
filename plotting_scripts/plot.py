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

    def plot_field_2d(self, plane='xy'):
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

            fx = self.bz
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

        # self.b = np.array([[[(x, y, z) for z in self.bz] for y in self.by] for x in self.bx])
    #
    # @staticmethod
    # def _find_nearest(array, value):
    #     array = np.asarray(array)
    #     idx = (np.abs(array - value).argmin())
    #     return idx
    #
    # def _get_plane_position(self, component, plane_position):
    #     if component == 'z':
    #         component = 2
    #         plane_idx = self._find_nearest(self.z_range, plane_position)
    #     elif component == 'y':
    #         component = 1
    #         plane_idx = self._find_nearest(self.y_range, plane_position)
    #     elif component == 'x':
    #         component = 0
    #         plane_idx = self._find_nearest(self.x_range, plane_position)
    #     else:
    #         component = 'abs'
    #         plane_idx = self._find_nearest(self.y_range, plane_position)
    #     return component, plane_idx
    #
    # @staticmethod
    # def _get_b_field_values_from_plane(b, component, plane_idx):
    #     if component == 0:
    #         return b[plane_idx]
    #     elif component == 1:
    #         return b[:, plane_idx]
    #     if component == 2:
    #         return b[:, :, plane_idx]
    #
    # def _get_b_field_values(self):
    #     return self.b
    #
    # def get_magnetic_field_value(self, component, plane_position):
    #     """Return the magnetic field value at a given plane."""
    #
    #     component, plane_idx = self._get_plane_position(component, plane_position)
    #
    #     if component == 'abs':
    #         b = self._get_b_field_values()
    #     else:
    #         b = self._get_b_field_values()[component]
    #
    #     return self._get_b_field_values_from_plane(b, component, plane_idx)
    #
    # def get_b(self, arg):
    #     x, y, z = arg
    #     local_b = self.b[(x, y, z)]
    #     return local_b
    #
    # def get_b_abs(self, arg):
    #     """
    #
    #     Parameters
    #     ----------
    #     arg: tuple
    #         x, y, z coordinates
    #     """
    #     # b = np.zeros(3)
    #     # print(self.b)
    #     local_b = self.b[arg]
    #     return np.linalg.norm(local_b)
    #
    # def get_2d_abs_plot_data(self, plane, source='new'):
    #     if source == 'storage':
    #         b, extent = read_data_from_file()
    #     else:
    #         b = self.get_magnetic_field_value(plane, plane_position=0)
    #
    #         if plane == 'yz':
    #             extent = (self.y_range.min(), self.y_range.max(), self.z_range.min(), self.z_range.max())
    #         elif plane == 'xy':
    #             extent = (self.x_range.min(), self.x_range.max(), self.y_range.min(), self.y_range.max())
    #         elif plane == 'xz':
    #             extent = (self.x_range.min(), self.x_range.max(), self.z_range.min(), self.z_range.max())
    #         else:
    #             extent = None
    #     return b, extent
    #
    # # def get_1d_b_values(self, component):
    # #     return self.get_magnetic_field_value(component='abs', plane_position=0)
    #
    # @staticmethod
    # def get_numerical_component(component):
    #     if component == 'x':
    #         return 0
    #     elif component == 'y':
    #         return 1
    #     elif component == 'z':
    #         return 2

    # def plot_field_1d_scalar(self, component):
    #     component = self.get_numerical_component(component)
    #
    #     plot_x_values = list()
    #     plot_y_values = list()
    #
    #     for point, b_field in self.b.items():
    #         x_value = point[0]
    #         if x_value not in plot_x_values:
    #             plot_x_values.append(x_value)
    #
    #         y_value = b_field[component]
    #         if y_value not in plot_y_values:
    #             plot_y_values.append(y_value)
    #
    #     fig = plt.figure()
    #     ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])  # main axes
    #
    #     ax.plot(plot_x_values, plot_y_values)
    #
    #     if self.x_ticks and self.x_ticks_labels:
    #         ax.set_xticks(self.x_ticks)
    #         ax.set_xticklabels(self.x_ticks_labels)
    #
    #     plt.show()
    #
    # def plot_field_1d_vec(self, _type='3d'):
    #     # Redefine plot ranges
    #     # y = len(self.x_range) * [0]
    #     # z = len(self.x_range) * [0]
    #     #
    #     # bx, by, bz = self.get_b_vec()
    #     # print(f'bx:{bx}\nby:{by}\nbz_{bz}')
    #     # x, y, z = np.meshgrid(self.x_range, self.y_range, self.z_range)
    #
    #     if _type == '3d':
    #         fig = plt.figure()
    #
    #         ax = fig.add_subplot(111, projection='3d')
    #
    #         for point, b_field in self.b.items():
    #             ax.quiver(point[0], point[1], point[2],
    #                       b_field[0], b_field[1], b_field[2],
    #                       # pivot='tip', length=vlength, arrow_length_ratio=0.3/vlength
    #                       )
    #
    #         ax.set_xlim([min(self.x_range), max(self.x_range)/3])
    #         ax.set_ylim([-0.01, 0.01])
    #         ax.set_zlim([-0.0025, 0.0025])
    #
    #         ax.set_xlabel('x')
    #         ax.set_ylabel('y')
    #         ax.set_zlabel('z')
    #
    #     elif _type == '2d':
    #         xi, yi = np.meshgrid(self.x_range, 0, indexing='ij')
    #         for x in self.x_range:
    #             v = self.b[x, 0, 0]
    #             print(v)
    #             # plt.axes([0.065, 0.065, 0.9, 0.9])
    #             plt.quiver(xi, yi, v[0], v[2], alpha=.5)
    #             plt.quiver(xi, yi, v[0], v[2], edgecolor='k', facecolor='none', linewidth=.5)
    #
    #     plt.show()
    #
    # def plot_field_2d_abs(self, plane=None):
    #     b, extent = self.get_2d_abs_plot_data(plane)
    #
    #     plt.imshow(b, aspect='auto', cmap=cm.magma, extent=extent)
    #
    #     plt.colorbar()
    #     plt.show()
    #
    # def get_2d_vec_plot_data(self):
    #     by = np.array([[[self.b[(x, y, z)][0] for z in self.z_range]
    #                    for y in self.y_range] for x in self.x_range])[:, 0]
    #
    #     bz = np.array([[[self.b[(x, y, z)][1] for z in self.z_range]
    #                    for y in self.y_range] for x in self.x_range])[:, :, 0]
    #     return by, bz

    # def plot_field_2d_vec(self):
    #     by, bz = self.get_2d_vec_plot_data()
    #
    #     plt.quiver(self.y_range, self.z_range, by, bz)
    #     plt.xlim(-1, 1)
    #     plt.ylim(-1, 1)
    #
    #     plt.show()


if __name__ == "__main__":
    plotter = Plotter()

    # plotter.plot_field_1d_scalar()
    plotter.plot_field_2d(plane='xy')
    # plotter.plot_field_3d()
