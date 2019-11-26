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

from utils.helper_functions import read_data_from_file, find_list_length_of_different_items
from experiments.mieze.parameters import absolute_x_position, HelmholtzSpinFlipper_position_HSF1, step, lambda_n


class Plotter:

    def __init__(self):
        self.x_range, self.y_range, self.z_range, self.bx, self.by, self.bz = read_data_from_file()

        self.preadjust_values()

    def preadjust_values(self):
        self.bx = 1e-1 * np.asarray(np.abs(self.bx))
        self.by = np.abs(self.by)

    def plot_field_1d_scalar(self, component=None):

        fig = plt.figure()
        # ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])  # main axes
        if component == 'x':
            plt.plot(self.x_range, self.bx)
        elif component == 'y':
            plt.plot(self.x_range, self.by)
        elif component == 'z':
            plt.plot(self.x_range, self.bz)
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

            f = self.bz
            # fy = self.by
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

    def extra_coil_check(self, b_extra, b_polariser):
        if b_extra.units != b_polariser.units:
            raise RuntimeError(f'B field units do not match:\n extra: {b_extra.units}\n'
                               f'polariser: {b_polariser.units}')

        if b_extra > b_polariser:
            pass
            # raise RuntimeError("Field of the extra coils too strong. Try resetting the extra coils.")

    def plot_adiabatic_check(self):
        index_first_hsf = np.argmin(abs(np.asarray(self.x_range)-HelmholtzSpinFlipper_position_HSF1))

        fig1, ax1 = plt.subplots()

        color = 'tab:red'
        ax1.set_xlabel('Neutron Trajectory (m)')
        ax1.set_ylabel('Magnetic field (G)', color=color)

        # logger.error(f'{absolute_x_position}\n]{Bx_values}\n{By_values}')

        ax1.plot(self.x_range, self.bx)
        ax1.plot(self.x_range, self.bz)

        ax1.legend((r'$B_x$', r'$B_z$'), loc=9)
        ax1.tick_params(axis='y', labelcolor=color)

        ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

        color = 'tab:green'
        x_pos = self.x_range[index_first_hsf]
        # noinspection PyTypeChecker
        theta_values = np.degrees(np.arctan(np.divide(np.asarray(self.bz), np.asarray(self.bx))))
        y_pos = theta_values[index_first_hsf]


        # we already handled the x-label with ax1
        ax2.set_ylabel(r'$\theta$ (degree)', color=color)
        ax2.plot(self.x_range, theta_values, color=color)
        ax2.tick_params(axis='y', labelcolor=color)
        ax2.plot(x_pos, y_pos, color='black', marker='o')
        ax2.text(x_pos, y_pos*0.9, '{:.1f}°'.format(y_pos))

        # fig1.tight_layout()  # otherwise the right y-label is slightly clipped
        plt.show()
        # plt.savefig('By_Bx.pdf')
        # plt.close()

        fig1, ax = plt.subplots()

        dtheta_dy = np.abs(np.gradient(theta_values, step))
        b_values = np.sqrt(np.power(np.asarray(self.bx), 2) + np.power(np.asarray(self.bz), 2))

        # ax.set_yscale('log')
        ax.plot(self.x_range, np.asarray(dtheta_dy) * 1e-2)  # y from m to cm
        ax.plot(self.x_range, 2.65*lambda_n * np.asarray(b_values) * 1e-1)  # B from Gauss to mT,
        ax.legend((r'$\frac{d\theta}{dy}$', r'$2.65\lambda B$'))

        ax.set_xlabel('Neutron Trajectory (m)')
        ax.set_ylabel("(degrees/cm)")
        ax.grid()
        # fig1.tight_layout()
        plt.show()
        # plt.savefig('Adiabatic_Check.pdf')
        # plt.close()


if __name__ == "__main__":
    plotter = Plotter()

    plotter.plot_field_1d_scalar(component='x')

    # plotter.plot_field_2d_clr_map(plane='yz')
    # plotter.plot_field_2d_vec_map(plane='xy')
    # plotter.plot_field_3d()

    # plotter.plot_adiabatic_check()
