# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019, 2020 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Generate an ideal magnetic field."""

import numpy as np

from simulation.parameters_simulation import absolute_x_position
from utils.helper_functions import save_obj, save_data_to_file, rotate


def main():
    """Compute and save the magnetic field."""

    data = dict()
    initial_b_field = np.array([0, 1, 0])
    for pos_x in absolute_x_position:
        angle = pos_x
        data[pos_x, 0, 0] = rotate(initial_b_field, angle, np.array([0, 0, 1]))

    save_data_to_file(data, '../data/ideal_magnetic_field_data.csv')
    save_obj(data, '../data/ideal_magnetic_field_data')


if __name__ == '__main__':
    main()
