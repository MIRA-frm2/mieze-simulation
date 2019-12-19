# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Plotting the magnetic field to check whether the adiabatic condition is fulfilled.

"""

import matplotlib.pyplot as plt
import numpy as np

from experiments.mieze.parameters import HelmholtzSpinFlipper_position_HSF1, lambda_n, step_x

from utils.helper_functions import read_data_from_file, append_column_to_csv
from utils.physics_constants import factor_T_to_G


def adjust_angle_gradient_for_adiabatic_condition(dtheta_dy):
    """Convert from angle/m to angle/cm"""
    return np.asarray(dtheta_dy) / 1e2


def compute_adiabatic_condition(b_values):
    """Compute the adiabatic condition value from the magnetic field.

    Note: 1 T = 1000 Gauss, hence 1 mT = 10 G

    Parameters
    ----------
    b_values: np.array, float
        Either the array of magnetic field values, or a single  value,

    Returns
    -------
    out: np.array, float
        The adiabatic transition condition value.
        Depending on the input, either the array or a single value.
    """
    factor_gauss_to_militesla = 1e3 / factor_T_to_G
    return 2.65 * lambda_n * np.asarray(b_values) * factor_gauss_to_militesla


def get_b_field_magnitude(bx, by):
    return np.sqrt(np.power(np.asarray(bx), 2) + np.power(np.asarray(by), 2))


class MyPlotter:
    """Customized Plotter Class."""

    def __init__(self, data_file='./data_adiabatic_transition.csv'):
        self.data_file = data_file
        self.x_range, self.y_range, self.z_range, self.bx, self.by, self.bz = read_data_from_file(data_file)

        self.theta_values = None
        self.dtheta_dy = None
        self.b_values = None
        self.adiabatic_condition_values = None

        self.preadjust_values()

    def preadjust_values(self):
        """Adjust values format."""
        self.bx = np.asarray(np.abs(self.bx))
        self.by = np.abs(self.by)

    def get_adiabatic_values(self, new=True):
        """Compute/Get the adiabatic specific values."""
        if new:
            self.theta_values = adjust_angle_gradient_for_adiabatic_condition(
                np.degrees(np.arctan(np.divide(np.asarray(self.by), np.asarray(self.bx)))))
            self.dtheta_dy = np.abs(np.gradient(self.theta_values, step_x))

            self.b_values = get_b_field_magnitude(self.bx, self.by)
            self.adiabatic_condition_values = compute_adiabatic_condition(self.b_values)

            append_column_to_csv(self.data_file, 'd theta/dy', self.dtheta_dy)
            append_column_to_csv(self.data_file, 'B magnitude', self.b_values)

    def plot_adiabatic_check(self):
        """Perform the adiabatic checks plots."""
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
        y_pos = self.theta_values[index_first_hsf]

        # we already handled the x-label with ax1
        ax2.set_ylabel(r'$\theta$ (degree)', color=color)
        ax2.plot(self.x_range, self.theta_values, color=color)
        ax2.tick_params(axis='y', labelcolor=color)
        ax2.plot(x_pos, y_pos, color='black', marker='o')
        ax2.text(x_pos, y_pos * 0.9, '{:.1f}°'.format(y_pos))

        # fig1.tight_layout()  # otherwise the right y-label is slightly clipped
        plt.show()
        # plt.savefig('By_Bx.pdf')
        # plt.close()

        fig1, ax = plt.subplots()

        ax.set_yscale('log')

        ax.plot(self.x_range, self.dtheta_dy)  # y from m to cm
        ax.plot(self.x_range, self.b_values)  # B from Gauss to mT,
        ax.legend((r'$\frac{d\theta}{dy}$', r'$2.65\lambda B$'))

        ax.set_xlabel('Neutron Trajectory (m)')
        ax.set_ylabel("(degrees/cm)")
        ax.grid()
        # fig1.tight_layout()
        plt.show()
        # plt.savefig('Adiabatic_Check.pdf')
        # plt.close()

    def get_adiabatic_difference_data(self):
        self.diff_data = self.adiabatic_condition_values - self.dtheta_dy

    def plot_adiabatic_check_difference(self):
        self.get_adiabatic_difference_data()
        plt.plot(self.x_range, self.diff_data)
        plt.yscale('log')
        plt.show()


if __name__ == "__main__":
    plot = MyPlotter()
    plot.get_adiabatic_values(new=True)
    # plot.plot_adiabatic_check()
    plot.plot_adiabatic_check_difference()
