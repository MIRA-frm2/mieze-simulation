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

from simulation.experiments.mieze.parameters import step_x, neutron_speed

from analysises.adiabatic_check.scripts.adiabatic_check import get_b_field_magnitude

from utils.helper_functions import read_data_from_file, convert_between_m_and_cm


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
    gyromagnetic_ratio = - 1.8324717143e8  # [rad/(s *T)]
    conversion_factor_angstrom_to_m = 1e10
    conversion_factor_cm_to_m = 1e2

    amplification_factor = 100

    prefactor = - gyromagnetic_ratio / (neutron_speed * conversion_factor_angstrom_to_m / conversion_factor_cm_to_m)
    prefactor *= amplification_factor
    return prefactor * b_field / d_theta


def main():
    """Plot the adiabatic parameter E."""
    data_file_magnetic_field = '../data/ideal_magnetic_field_data.csv'
    x_range_bfield, y_range, z_range, bx, by, bz = read_data_from_file(data_file_magnetic_field)
    data_file_polarisation = '../data/data_polarisation.csv'
    x_range_pol, y_range, z_range, pol_x, pol_y, pol_z = read_data_from_file(data_file_polarisation)

    dtheta_dy_values = compute_dtheta_dy(bx, by)
    b_values = get_b_field_magnitude(bx, by)
    adiabatic_parameter_e = list(compute_adiabatic_parameter_e(b_values, dtheta_dy_values))

    polarisation_data = list()

    n = len(x_range_pol)
    for i in range(n):
        i = n - 1 - i
        b_vec = np.array([bx[i], by[i], bz[i]])
        pol_vec = np.array([pol_x[i], pol_y[i], pol_z[i]])

        normalisation = np.linalg.norm(pol_vec) * np.linalg.norm(b_vec)
        print(f'b_vec: {b_vec}\npol_vec: {pol_vec}')
        polarisation = np.dot(b_vec, pol_vec) / normalisation

        polarisation_data.append(polarisation)

        if x_range_bfield[i] not in x_range_pol:
            adiabatic_parameter_e.pop(i)
    adiabatic_parameter_e.pop(0)

    fig, ax1 = plt.subplots()

    color = 'tab:red'
    ax1.set_xlabel('x beamline [m]')
    ax1.set_ylabel('Polarisation', color=color)
    ax1.plot(x_range_pol, polarisation_data, color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:blue'
    ax2.set_ylabel('Adiabatic Parameter E', color=color)  # we already handled the x-label with ax1
    ax2.plot(x_range_pol, adiabatic_parameter_e, color=color)
    ax2.tick_params(axis='y', labelcolor=color)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()


if __name__ == '__main__':
    main()
