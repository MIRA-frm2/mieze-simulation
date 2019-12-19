# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Helper functions used for various purposes."""

import csv
import logging
import numpy as np
import pandas as pd
import pickle

# Create a custom logger
logger = logging.getLogger(__name__)


def adjust_field(vector):
    """

    Parameters
    ----------
    vector: np.array

    """
    return np.array([vector[2], vector[1], vector[0]])


def append_column_to_csv(filename, column_name, column_data):
    csv_input = pd.read_csv(filename)
    csv_input[column_name] = column_data
    csv_input.to_csv(filename, index=False)


def convert_between_m_and_cm(value, backwards=False):
    """Convert from m to cm, and backwards.

    Parameters
    ----------
    value: np.array, float, ndarray
        Array or value to be converted.
    backwards: bool
        Flag indicating whether to convert from m to cm, or from cm to m.
        If True, converts from m to cm.
        If False, converts from cm to m.
    """
    if backwards:
        return np.asarray(value) * 1e2
    else:
        return np.asarray(value) / 1e2


def find_list_length_of_different_items(x):
    """Returns the length if different items from a list.

    The initial list may contain duplicates that would otherwsid
    """
    xx = list()
    for item in x:
        if item not in xx:
            xx.append(item)
    return len(xx)


def find_nearest(array, value, index=True):
    """Find the nearest grid point index from the requested value.

    Parameters
    ----------
    array: np.array
        Array with values.
    value: float
        Float to be found close an array element.
    index: bool, optional
        Flag indicating whether to retrieve value as the index, or the value of the array at that index.

    Returns
    -------
    idx: int
        The index of the array element closest to the value.

    >>> find_nearest([3, 2, 1], 1.1)
    2
    >>> find_nearest([3, 2, 1], 1.1, index=False)
    1

    """
    array = np.asarray(array)
    idx = (np.abs(array - value).argmin())
    if index:
        return idx
    else:
        return array[idx]


def get_phi(y, z):
    """Compute the angle from two coordinates."""
    if y:
        return np.arctan(z/y)
    else:
        if z >= 0:
            return np.pi / 2
        else:
            return np.pi * 3 / 2


def load_obj(name):
    """Load object (magnetic field) as pickled object."""
    with open(f'{name}.pkl', 'rb') as f:
        return pickle.load(f)


def read_data_from_file(file_name):
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
                # print(f'Column names are {", ".join(row)}')
                line_count += 1
            else:
                try:
                    x.append(float(row[0]))
                    y.append(float(row[1]))
                    z.append(float(row[2]))

                    bx.append(float(row[3]))
                    by.append(float(row[4]))
                    bz.append(float(row[5]))
                except IndexError:
                    bx.append(float(row[1]))

        return x, y, z, bx, by, bz


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


def rotate(vector, phi, axis):
    """Rotate the vector with an angle phi with respect to the axis."""
    n = axis / np.linalg.norm(axis)
    c = np.cos(phi)
    s = np.sin(phi)

    n1 = n[0]
    n2 = n[1]
    n3 = n[2]

    r = [[n1 ** 2 * (1 - c) + c, n1 * n2 * (1 - c) - n3 * s, n1 * n3 * (1 - c) + n2 * s],
         [n2 * n1 * (1 - c) + n3 * s, n2 ** 2 * (1 - c) + c, n2 * n3 * (1 - c) - n1 * s],
         [n3 * n1 * (1 - c) - n2 * s, n3 * n2 * (1 - c) + n1 * s, n3 ** 2 * (1 - c) + c]]
    return np.dot(r, vector)


def save_data_to_file(data, file_name, extension='.csv'):
    """Save data to file."""
    if extension in file_name:
        full_filename = file_name
    else:
        full_filename = f'{file_name}{extension}'

    with open(full_filename, 'w') as file:

        logger.info(f'Writing data to file {full_filename}')

        csv_writer = csv.writer(file, delimiter=',')
        csv_writer.writerow(["x", "y", "z", "Bx", "By", "Bz"])

        for point, field in data.items():
            # logger.debug(type(point))
            if type(point) != np.float64:
                point = list(point)
                row = point + list((field[0], field[1], field[2]))
            else:
                row = list([point, field])
            csv_writer.writerow(row)


def save_obj(obj, name):
    """Save object (magnetic field) as pickled object."""
    with open(f'{name}.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


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
