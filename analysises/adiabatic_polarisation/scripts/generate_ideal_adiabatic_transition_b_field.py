# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Generate an ideal magnetic field."""

import numpy as np

from simulation.experiments.mieze.parameters import absolute_x_position
from utils.helper_functions import save_obj, save_data_to_file


def compute_ideal_b_field(x):
    """Compute an ideal magnetic field.

    Parameters
    ----------
    x: float, ndarray

    Returns
    -------
    out: float, ndarray
    """
    strong_b_field = 10
    return np.array([20 * x, strong_b_field - 20 * x, 0])


def main():
    """Compute and save the magnetic field."""

    data = dict()
    for pos_x in absolute_x_position:
        data[pos_x, 0, 0] = compute_ideal_b_field(pos_x)

    save_data_to_file(data, '../data/ideal_magnetic_field_data.csv')
    save_obj(data, '../data/ideal_magnetic_field_data')


if __name__ == '__main__':
    main()
