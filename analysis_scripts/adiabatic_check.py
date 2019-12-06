# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Plotting the magnetic field to check whether the adiabatic condition is fulfilled."""

import matplotlib.pyplot as plt
import numpy as np

from experiments.mieze.main import Mieze
from elements.coils import Coil
from experiments.mieze.parameters import HelmholtzSpinFlipper_position_HSF1, lambda_n, step
from experiments.mieze.parameters import L1, L2

from utils.helper_functions import read_data_from_file, save_data_to_file


class Plotter:

    def __init__(self):
        self.x_range, self.y_range, self.z_range, self.bx, self.by, self.bz = read_data_from_file()

        self.preadjust_values()

    def preadjust_values(self):
        self.bx = np.asarray(np.abs(self.bx))
        self.by = np.abs(self.by)

    def plot_adiabatic_check(self):
        index_first_hsf = np.argmin(abs(np.asarray(self.x_range) - HelmholtzSpinFlipper_position_HSF1))

        fig1, ax1 = plt.subplots()

        color = 'tab:red'
        ax1.set_xlabel('Neutron Trajectory (m)')
        ax1.set_ylabel('Magnetic field (G)', color=color)

        # logger.error(f'{absolute_x_position}\n]{Bx_values}\n{By_values}')

        ax1.plot(self.x_range, self.bx)
        ax1.plot(self.x_range, self.by)

        ax1.legend((r'$B_x$', r'$B_y$'), loc=9)
        ax1.tick_params(axis='y', labelcolor=color)

        ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

        color = 'tab:green'
        x_pos = self.x_range[index_first_hsf]
        # noinspection PyTypeChecker
        theta_values = np.degrees(np.arctan(np.divide(np.asarray(self.by), np.asarray(self.bx))))
        y_pos = theta_values[index_first_hsf]

        # we already handled the x-label with ax1
        ax2.set_ylabel(r'$\theta$ (degree)', color=color)
        ax2.plot(self.x_range, theta_values, color=color)
        ax2.tick_params(axis='y', labelcolor=color)
        ax2.plot(x_pos, y_pos, color='black', marker='o')
        ax2.text(x_pos, y_pos * 0.9, '{:.1f}°'.format(y_pos))

        # fig1.tight_layout()  # otherwise the right y-label is slightly clipped
        plt.show()
        # plt.savefig('By_Bx.pdf')
        # plt.close()

        fig1, ax = plt.subplots()

        dtheta_dy = np.abs(np.gradient(theta_values, step))
        b_values = np.sqrt(np.power(np.asarray(self.bx), 2) + np.power(np.asarray(self.by), 2))

        # ax.set_yscale('log')
        ax.plot(self.x_range, np.asarray(dtheta_dy) * 1e-2)  # y from m to cm
        ax.plot(self.x_range, 2.65 * lambda_n * np.asarray(b_values) * 1e-1)  # B from Gauss to mT,
        ax.legend((r'$\frac{d\theta}{dy}$', r'$2.65\lambda B$'))

        ax.set_xlabel('Neutron Trajectory (m)')
        ax.set_ylabel("(degrees/cm)")
        ax.grid()
        # fig1.tight_layout()
        plt.show()
        # plt.savefig('Adiabatic_Check.pdf')
        # plt.close()


if __name__ == "__main__":
    plot = Plotter()
    plot.plot_adiabatic_check()
