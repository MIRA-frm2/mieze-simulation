# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019, 2020 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""An analysis of the polarizer to better fit the simulation to the physical instrument."""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from simulation.elements.polariser import Polariser

# These values represent measured data
# The 0 actually corresponds to the right side of the polariser
length = [0, 0.025, 0.05, 0.075, 0.10, 0.20, 0.30]  # [m]
b = np.asarray([52.5, 24, 9, 3.7, 2.1, 0.4, 0.12]) * 10  # [G]


def fit_polariser_data(polariser):
    """Perform the fit for the polariser data."""
    position = np.array(length)
    b_field_values = np.array(b)

    popt, pcov = curve_fit(polariser.b_field_fitted, xdata=position, ydata=b_field_values, p0=[3, 700])
    return popt, pcov


def check():
    polariser = Polariser(position=0)

    popt, pcov = fit_polariser_data(polariser)

    computed_values = list()
    for pos in length:
        computed_values.append(np.abs(polariser.b_field([pos, 0, 0])))

    fig, ax = plt.subplots()
    plt.plot(length, computed_values)
    plt.plot(length, b, '*')
    plt.plot(length, polariser.b_field_fitted(length, power=popt[0], amplitude=popt[1]))

    ax.legend((None, r'Theoretical', None, r'Measured', r'Fitted'), loc=9)

    ax.set_xlabel('Neutron Trajectory (m)')
    ax.set_ylabel('Magnetic field (G)')

    plt.yscale('log')
    plt.show()


if __name__ == "__main__":
    check()
