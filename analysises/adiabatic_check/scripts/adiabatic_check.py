# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019, 2020 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Plotting the magnetic field to check whether the adiabatic condition is fulfilled."""

import matplotlib.pyplot as plt
import numpy as np

from experiments.mieze.parameters import HelmholtzSpinFlipper_position_HSF1
from simulation.parameters_simulation import step_x
from simulation.beamline.beamline_properties import wavelength

from utils.helper_functions import append_column_to_csv, convert_between_m_and_cm, read_data_from_file
from utils.physics_constants import factor_T_to_G


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
    return 2.65 * wavelength * np.asarray(b_values) * factor_gauss_to_militesla


def get_b_field_magnitude(bx, by, bz=None):
    """Compute the magnetic field magnitude from x and y components only.

    Parameters
    ----------
    bx: ndarray, float
        Magnetic field on x axis.
    by: ndarray, float
        Magnetic field on y axis.
    bz: ndarray, float, optional
        Magnetic field on z axis.
        By default, set to None, to compute the magnitude from only two components.

    Returns
    -------
    out: np.array, float
        Magnetic field magnitude computed from x and y components only.
    """
    if bz is not None:
        return np.sqrt(np.power(np.asarray(bx), 2) + np.power(np.asarray(by), 2) + np.power(np.asarray(bz), 2))
    else:
        return np.sqrt(np.power(np.asarray(bx), 2) + np.power(np.asarray(by), 2))


class MyPlotter:
    """Customized Plotter Class."""

    def __init__(self, data_file='../../data/data_magnetic_field.csv'):
        """

        Parameters
        ----------
        data_file: string
            Location of the data file containing the magnteic field values.
        """
        self.data_file = data_file
        self.x_range, self.y_range, self.z_range, self.bx, self.by, self.bz = read_data_from_file(data_file)

        self.theta_values = None
        self.dtheta_dy = None
        self.b_values = None
        self.adiabatic_condition_values = None
        self.diff_data = None

        self.preadjust_values()

    def preadjust_values(self):
        """Adjust values format for plotting/analysis purposes."""
        self.bx = np.asarray(np.abs(self.bx))
        self.by = np.abs(self.by)

    def get_adiabatic_values(self, new=True):
        """Compute/Get the adiabatic specific values."""
        if new:
            self.theta_values = convert_between_m_and_cm(
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

        ax1.set_title(f'Helmholtz Coils at position: {round(HelmholtzSpinFlipper_position_HSF1 *100, 2)} cm')

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

        ax.plot(self.x_range, self.dtheta_dy)
        ax.plot(self.x_range, self.b_values)

        ax.set_title(f'Helmholtz Coils at position: {round(HelmholtzSpinFlipper_position_HSF1 *100, 2)} cm')
        ax.legend((r'$\frac{d\theta}{dy}$', r'$2.65\lambda B$'))

        ax.set_xlabel('Neutron Trajectory (m)')
        ax.set_ylabel("(degrees/cm)")
        ax.grid()
        plt.show()
        # plt.savefig('Adiabatic_Check.pdf')
        # plt.close()

    def get_adiabatic_difference_data(self):
        """Compute the adiabatic transition discriminant as the difference between the two values."""
        self.diff_data = self.adiabatic_condition_values - self.dtheta_dy

    def plot_adiabatic_check_difference(self):
        """Plot the adiabatic transition discriminant."""
        self.get_adiabatic_difference_data()

        plt.plot(self.x_range, self.diff_data)
        plt.yscale('log')
        plt.show()
