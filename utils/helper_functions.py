# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Helper functions used for various purposes."""

import numpy as np
import csv


def transform_cartesian_to_cylindrical(x, y, z):
    """Transform coordinates from cartesian to cylindrical."""
    x = x
    rho = np.sqrt(y ** 2 + z ** 2)
    if z:
        theta = np.arctan(y / z)
    else:
        theta = 0
    return x, rho, theta


def transform_cylindrical_to_cartesian(x, rho, theta):
    """Transform coordinates from cylindrical to cartesian."""
    x = x
    y = rho * np.cos(theta)
    z = rho * np.sin(theta)
    return x, y, z


def get_phi(y, z):
    """Compute the angle from two coordinates."""
    if y:
        return np.arctan(z/y)
    else:
        if z >= 0:
            return np.pi / 2
        else:
            return np.pi * 3 / 2


def adjust_field(vector):
    """

    Parameters
    ----------
    vector: np.arry

    """
    return np.array([vector[2], vector[1], vector[0]])


def save_data_to_file(data, file_name, extension='.csv'):
    """Save data to file."""
    full_filename = f'{file_name}{extension}'

    with open(full_filename, 'w') as file:

        print(f'Writing data to file {full_filename}')

        csv_writer = csv.writer(file, delimiter=',')
        csv_writer.writerow(["x", "y", "z", "Bx", "By", "Bz"])

        for point, field in data.items():
            row = list(point) + list(field)
            csv_writer.writerow(row)


def read_data_from_file():
    pass
