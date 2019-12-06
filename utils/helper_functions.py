# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Helper functions used for various purposes."""

import logging
import numpy as np
import csv
import pickle

# Create a custom logger
logger = logging.getLogger(__name__)


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
    vector: np.array

    """
    return np.array([vector[2], vector[1], vector[0]])


def save_obj(obj, name):
    """Save object (magnetic field) as pickled object."""
    with open(f'{name}.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def load_obj(name):
    """Load object (magnetic field) as pickled object."""
    with open(f'{name}.pkl', 'rb') as f:
        return pickle.load(f)


def save_data_to_file(data, file_name, extension='.csv'):
    """Save data to file."""
    full_filename = f'{file_name}{extension}'

    with open(full_filename, 'w') as file:

        logger.info(f'Writing data to file {full_filename}')

        csv_writer = csv.writer(file, delimiter=',')
        csv_writer.writerow(["x", "y", "z", "Bx", "By", "Bz"])

        # print(data)
        for point, field in data.items():
            # logger.debug(type(point))
            if type(point) != np.float64:
                point = list(point)
            else:
                point = list((point, ))

            # print(point, field)

            row = point + list((field[0], field[1], field[2]))
            csv_writer.writerow(row)


def read_data_from_file(file_name='../data/data.csv'):
    """Read data from file"""
    x = list()
    y = list()
    z = list()

    bx = list()
    by = list()
    bz = list()

    with open(file_name) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                print(f'Column names are {", ".join(row)}')
                line_count += 1
            else:
                x.append(float(row[0]))
                y.append(float(row[1]))
                z.append(float(row[2]))

                bx.append(float(row[3]))
                by.append(float(row[4]))
                bz.append(float(row[5]))

    return x, y, z, bx, by, bz


def find_list_length_of_different_items(x):
    """Returns the length if different items from a list.

    The initial list may contain duplicates that would otherwsid
    """
    xx = list()
    for item in x:
        if item not in xx:
            xx.append(item)
    return len(xx)


def sanitize_output(func):
    def wrapper_sanitize_output(*args, **kwargs):
        value = func(*args, **kwargs)
        if abs(value) > 10e4:
            value = 0
        return value
    return wrapper_sanitize_output


def unit_square(x_min, x_max, grid):
    """Return an array of a unit square function.

    """
    x_values = list()
    for x in grid:
        if x_min < x < x_max:
            x_values.append(1)
        else:
            x_values.append(0)
    return np.array(x_values)
