# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019, 2020 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Over layer that performs the analysis."""

import numpy as np

from analysises.adiabatic_check.scripts.find_relative_distance import find_relative_distance
from analysises.adiabatic_check.scripts.adiabatic_check import MyPlotter


def find_ideal_relative_distance():
    iteration_values = np.arange(0.0, 0.2 + 0.025, step=0.025)

    find_relative_distance(iteration_values)


def plot_adiabatic_condition():
    plot = MyPlotter()
    plot.get_adiabatic_values(new=True)
    plot.plot_adiabatic_check()


if __name__ == "__main__":
    # find_ideal_relative_distance()
    plot_adiabatic_condition()
