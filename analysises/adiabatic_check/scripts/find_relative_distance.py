# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Search for the ideal relative distance between polariser and helmholtz coils."""

import matplotlib.pyplot as plt

from simulation.experiments.mieze.main_mieze import main_mieze
from analysises.adiabatic_check.scripts.adiabatic_check import MyPlotter

from utils.helper_functions import convert_between_m_and_cm


def find_relative_distance(iteration_values):

    diff_data_list = list()

    # Compute the adiabatic transition values for each iteration
    for val in iteration_values:
        data_file = f'./data/test_data_{int(val*100)}.csv'

        main_mieze(filename=data_file, spin_flipper_distance=val)

        # Using the existing Plotter class, read values and compute values
        plottter = MyPlotter(data_file=data_file)
        plottter.get_adiabatic_values(new=True)
        plottter.get_adiabatic_difference_data()

        diff_data_list.append(plottter.diff_data)

    # Plot the values
    for item in diff_data_list:
        plt.plot(plottter.x_range, item)

    # Plot customization
    plt.xlabel('Neutron Trajectory (m)')
    plt.ylabel(r'Adiabatic condition discriminant: $2.65 \cdot B \cdot \lambda - \dfrac{d\theta}{dy}$')

    plt.yscale('log')
    plt.legend(convert_between_m_and_cm(iteration_values, backwards=True),  # Convert units from [m] to [cm]
               title="Relative distance between polariser and helmholtz coils [cm]",
               loc=1)  # Place the legend in the upper right corner
    plt.show()
