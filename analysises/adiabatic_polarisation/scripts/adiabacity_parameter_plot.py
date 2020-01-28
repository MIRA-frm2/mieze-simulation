# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Plot the adiabatic parameter E."""

import matplotlib.pyplot as plt
import numpy as np

from simulation.experiments.mieze.parameters import neutron_speed, absolute_x_position, step_x
from simulation.particles.neutron import Neutron

from analysises.adiabatic_check.scripts.adiabatic_check import get_b_field_magnitude

from utils.helper_functions import convert_between_m_and_cm, rotate


def compute_dtheta_dy(bx, by):
    """Compute the angle gradient.

    Parameters
    ----------
    bx: float, ndarray
    by: float, ndarray

    Returns
    -------
    dtheta_dy: float, ndarray
    """
    theta_values = convert_between_m_and_cm(
        np.degrees(np.arctan(np.divide(np.asarray(by), np.asarray(bx)))))
    dtheta_dy = np.abs(np.gradient(theta_values, step_x))
    return dtheta_dy


def compute_adiabatic_parameter_e(b_field, d_theta):
    """Compute the adiabatic parameter E.

    Parameters
    ----------
    b_field: float, ndarray
    d_theta: float, ndarray

    Returns
    -------
    out: float, ndarray
    """
    gyromagnetic_ratio = - 1.8324717143e8  # [rad/(s * T)]
    conversion_factor_radian_to_degrees = 180 / np.pi
    conversion_factor_angstrom_to_m = 1e10
    conversion_factor_cm_to_m = 1e2

    amplification_factor = 200

    prefactor = - (gyromagnetic_ratio * conversion_factor_radian_to_degrees) \
        / (neutron_speed * conversion_factor_angstrom_to_m / conversion_factor_cm_to_m)
    prefactor *= amplification_factor
    return prefactor * b_field / d_theta


def compute_polarisation(polarisation_vector, magnetic_field_vector):
    """Compute the polarisation from the vectors.

    Parameters
    ----------
    polarisation_vector: ndarray
        3D polarisation vector
    magnetic_field_vector: ndarray
        3D magnetic field vector

    Returns
    -------
    polarisation: float
        Value between -1 and 1, specifying the polarisation.
    """
    normalisation = np.linalg.norm(polarisation_vector) * np.linalg.norm(magnetic_field_vector)
    polarisation = np.dot(magnetic_field_vector, polarisation_vector) / normalisation
    return polarisation


def plot_polarisation_and_adiabatic_parameter_against_beamline(polarisation_data, adiabatic_parameter_e):
    fig, ax1 = plt.subplots()

    color = 'tab:red'
    ax1.set_xlabel('x beamline [m]')
    ax1.set_ylabel('Polarisation', color=color)
    ax1.plot(absolute_x_position, polarisation_data, color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:blue'
    ax2.set_ylabel('Adiabatic Parameter E', color=color)  # we already handled the x-label with ax1
    ax2.plot(absolute_x_position, adiabatic_parameter_e, color=color)
    ax2.tick_params(axis='y', labelcolor=color)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()


def plot_polarisation_against_adiabatic_parameter(adiabatic_parameter_e_data, final_polarisation_data):
    plt.plot(adiabatic_parameter_e_data, final_polarisation_data, color='tab:red')
    plt.xlabel('Adiabatic parameter (E)')
    plt.ylabel('Neutron Polarisation')
    plt.show()


def main():
    """Plot the adiabatic parameter E."""
    final_polarisation_data = list()
    adiabatic_parameter_e_data = list()

    magnetic_field_values = np.linspace(0, 10, 200)

    for max_b_value in magnetic_field_values:
        initial_b_field = np.array([0, max_b_value, 0])
        bx, by, bz = list(), list(), list()
        initial_polarisation = np.array([0, 1, 0])

        polarisation_data = list()

        neutron = Neutron(velocity=np.array([neutron_speed, 0, 0]),
                          position=np.array([0, 0, 0]), polarisation=initial_polarisation)

        for pos_x in absolute_x_position:
            angle = pos_x * 10
            magnetic_field_vector = rotate(initial_b_field, angle, np.array([0, 0, 1]))
            # print(f'magnetic field vector: {magnetic_field_vector}')

            bx.append(magnetic_field_vector[0])
            by.append(magnetic_field_vector[1])
            bz.append(magnetic_field_vector[2])

            neutron.set_position_x(pos_x)
            cell_time = step_x/neutron_speed
            polarisation_vector = neutron.compute_polarisation(magnetic_field_vector=magnetic_field_vector,
                                                               time=cell_time)

            polarisation = compute_polarisation(polarisation_vector, np.asarray(magnetic_field_vector))

            polarisation_data.append(polarisation)

        dtheta_dy_values = compute_dtheta_dy(np.asarray(bx), np.asarray(by))
        b_values = get_b_field_magnitude(np.asarray(bx), np.asarray(by), np.asarray(bz))
        adiabatic_parameter_e = list(compute_adiabatic_parameter_e(b_values, dtheta_dy_values))
        # Repair the adiabatic parameter list as the first two values are badly computed
        adiabatic_parameter_e[0:2] = adiabatic_parameter_e[-3:-1]

        # Assign final values
        adiabatic_parameter_e_data.append(adiabatic_parameter_e[-1])
        final_polarisation_data.append(polarisation_data[-1])

        # print('Finished computing one magnetic field value.')

    # plot_polarisation_and_adiabatic_parameter_against_beamline(polarisation_data, adiabatic_parameter_e)
    plot_polarisation_against_adiabatic_parameter(adiabatic_parameter_e_data, final_polarisation_data)


if __name__ == '__main__':
    main()
