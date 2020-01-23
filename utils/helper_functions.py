# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Helper functions used for various purposes."""

import csv
import json
import logging
import numpy as np
import pandas as pd
import pickle


from utils.physics_constants import earth_field

# Create a custom logger
logger = logging.getLogger(__name__)


def adjust_field(vector):
    """

    Parameters
    ----------
    vector: np.array

    """
    return np.array([vector[2], vector[1], vector[0]])


def add_earth_magnetic_field(field, flag):
    """Add the earth magnetic field to the computation.

    Parameters
    ----------
    field: np.array, ndarray
        Magnetic field array.
    flag: bool
        Flag indicating whether to consider the earth magnetic field or not.
        If True, adds the earth magnetic field to the input field.
        If False, only returns the input field.

    >>> add_earth_magnetic_field(np.array([0, 0, 0]), True)
    array([ 0.   ,  0.21 , -0.436])
    >>> add_earth_magnetic_field(np.array([1, 2, 0]), True)
    array([ 1.   ,  2.21 , -0.436])
    >>> add_earth_magnetic_field(np.array([0, 0, 0]), False)
    array([0, 0, 0])
    >>> add_earth_magnetic_field(np.array([1, 2, 3]), False)
    array([1, 2, 3])

    """
    if flag:
        return field + earth_field
    else:
        return field


def append_column_to_csv(filename, column_name, column_data):
    """Append the column data with column name to the csv file filename.

    Parameters
    ----------
    filename: str
        Name of the file to append to.
    column_name: str
        Column name to append.
    column_data: str
        Column data to append.
    """
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

    >>> convert_between_m_and_cm(1, backwards=False)
    0.01
    >>> convert_between_m_and_cm(1, backwards=True)
    100.0

    """
    if backwards:
        return np.asarray(value) * 1e2
    else:
        return np.asarray(value) / 1e2


def find_list_length_of_different_items(x):
    """Returns the length of different items from a list.

    The initial list may contain duplicates that would otherwise be doubly counted for.

    >>> find_list_length_of_different_items([0, 1])
    2
    >>> find_list_length_of_different_items([0, 1, 3])
    3
    >>> find_list_length_of_different_items([0, 1, 1])
    2

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
    """Compute the angle from two coordinates.

    >>> get_phi(0, 1)
    1.5707963267948966
    >>> get_phi(0, -1)
    4.71238898038469
    >>> get_phi(1, 1)
    0.7853981633974483

    """
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
    """Read data from file.

    Parameters
    ----------
    file_name: str
        Name of the file to be read from.
    """
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
    """Define a wrapper to sanitie ouputs from infinty values."""
    def wrapper_sanitize_output(*args, **kwargs):
        value = func(*args, **kwargs)
        if abs(value) > 10e4:
            value = 0
        return value
    return wrapper_sanitize_output


def unit_square(x_min, x_max, grid):
    """Return an array of a unit square function.

    Parameters
    ----------
    x_min: int, float
        Left edge of the flat top/unit square.
    x_max: int, float
        Right edge of the flat top/unit square.
    grid: np.array
        Array of x values to be transformed to a unit square.

    Returns
    -------
    y_values: np.array
        Array of 0 or 1 corresponding to a flat top.

    >>> unit_square(0.1, 0.2, np.array([0.0, 0.1, 0.15, 0.2, 0.3]))
    array([0, 1, 1, 1, 0])

    """
    y_values = list()
    for x in grid:
        if x_min <= x <= x_max:
            y_values.append(1)
        else:
            y_values.append(0)
    return np.array(y_values)


def rotate(vector, phi, axis):
    """Rotate the vector with an angle phi with respect to the axis.

    Parameters
    ----------
    vector: ndarray, list
        Vector to be rotated.
    phi: float
        Angle to be rotate by, in radians, e.g. a rotation of np.pi/2 corresponds to 90 degrees.
    axis: np.array, list
        Axis to be rotated from.

    Returns
    -------
    out: ndarray
        Rotated vector by phi with respect to axis.

    >>> rotate([1, 0, 0], np.pi/2, [0, 1, 0])
    array([ 6.123234e-17,  0.000000e+00, -1.000000e+00])
    >>> rotate([1, 0, 0], np.pi, [0, 1, 0])
    array([-1.0000000e+00,  0.0000000e+00, -1.2246468e-16])
    """
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


def save_metadata_to_file(filename, metadata):
    """Write metadata to file.

    Parameters
    ----------
    filename: str
        Filename for data storage.
    metadata: dict
        Dictionary containing the metadata to be written.
    """
    with open(f'{filename}_metadata.txt', 'w') as f:
        description = f'The {filename} file has been generated using the following parameters: \n'
        f.write(description)
        for key, val in metadata.items():
            f.write(f'\n{key}\n')
            f.write(json.dumps(val))


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
