# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019, 2020 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Investigation script to find the ideal position for the two Coil Set Pairs."""

import matplotlib.pyplot as plt
from multiprocessing import Pool
import numpy as np
from scipy.stats import chisquare

from simulation.elements.coil_set import CoilSet
from experiments.mieze.parameters import LENGTH_COIL_INNER, DISTANCE_BETWEEN_INNER_COILS


from utils.helper_functions import unit_square, save_data_to_file


def minimizer_function(computed_values, expected_values, use_chisquare=False):
    """Evaluate the similarity of the computed and expected arrays.

    Parameters
    ----------
    computed_values: np.array
        The array of computed values.
    expected_values: np.array
        The array of ideal/expeted values, used as reference.
    use_chisquare: bool, optional
        Flag indicating whether to use the chisquare test.
        In order to use chisquare, one needs arrays with at least a value of 5 in each array element.
        Defaults to False.

    Returns
    -------
    fit_value: float
        Value indicating the similarity of the computed array to the expected one
    """
    if use_chisquare:
        fit_value = chisquare(1e4 * np.asarray(computed_values),
                              1e4 * np.asarray(expected_values))
    else:
        fit_value = 0
        n = len(computed_values)
        for j in range(n):
            dif = computed_values[j] - expected_values[j]
            fit_value += dif ** 2
        fit_value = 1 / (fit_value * n)

    return fit_value


def get_ideal_field(b_max, middle_position, x_positions, scale_factor=0.99):
    """Returns the ideal magnetic field.

    Parameters
    ----------
    b_max: float
        Flat TOP value for the magnetic field.
        It is usually the maximum of the gaussian peak
    middle_position: float
        Position  of the centre of the ideal field.
        Should be the centre in between the Helmholtz coils.
    x_positions: np.array
        Array of positions where the ideal field is supposed to be computed at.
    scale_factor: float, optional
        Perhaps the ideal magnetic field is not supposed to be at the maximum of the gaussian peak.
        This factor allows to implement a scale that reduces it.
        Should have values cloes to unity.
        Defaults to 0.99.

    Returns
    -------
    out: np.array
        Array representing the ideal magnetic field value.
    """
    return scale_factor * b_max \
           * unit_square(middle_position - DISTANCE_BETWEEN_INNER_COILS/2 - LENGTH_COIL_INNER - LENGTH_COIL_INNER * 0.5,
                         middle_position + DISTANCE_BETWEEN_INNER_COILS / 2 + LENGTH_COIL_INNER * 0.5,
                         x_positions)


def define_iteration_values_for_coil_distances(n=1, max_distance=0.01):
    """Define the iteration values.

    Parameters
    ----------
    n: int, optional
        Number of iteration steps.
        Defaults to 1.
    max_distance: float, optional
        The maximum distance to be iterated over.
        Defaults to 0.01.

    Returns
    -------
    out: np.array
        Array consisting of iteration points for the analysis.
    """
    if n == 1:
        lspace = (0.0,)
    else:
        lspace = np.linspace(0, max_distance, n)

    return lspace


def define_computational_grid():
    """Define the computational grid space."""
    start_point = -0.35  # [m]
    end_point = 0.35  # [m]  # Positions.get_position_coilA()
    return np.linspace(start_point, end_point, num=700)


def plot_l1l2_cmap(fits, l, plot_name=None):
    """Plot the 2d color map."""
    if len(l) > 1:
        plt.imshow(fits, aspect='auto', extent=[min(l), max(l), min(l), max(l)], origin='lower')
        plt.colorbar()
        plt.xlabel('L1 (First Outer Coil Distance) [m]')
        plt.ylabel('L2 (Second Outer Coil Distance) [m]')

        # plt.show()
        plt.savefig(f'../docs/experiments/MIEZE/coil_set_optimization{plot_name}.png')
        plt.close()


def plot_ideal_position(middle_position, best_fit, x_positions, coil_type=None, plot_name=None):
    """Plot the magnetic field for the ideal position."""

    # Plot values
    print(best_fit)
    coil_set = CoilSet(coil_type=coil_type,
                       name='CoilSet',
                       position=middle_position,
                       distance_12=best_fit['l12'],
                       distance_34=best_fit['l34'])

    # Compute B_field_values
    with Pool(4) as p:
        b_field_values = p.map(coil_set.b_field, x_positions)

    dictionary = dict(zip(x_positions, b_field_values))
    save_data_to_file(dictionary, file_name='ideal_position_b_field_values')

    fig = plt.figure()
    ax = plt.axes()

    plt.plot(x_positions, b_field_values)

    # Plot idea B field
    b_max = max(b_field_values)
    ideal_b_field = get_ideal_field(b_max, middle_position, x_positions)

    # plt.plot(x_positions, ideal_b_field)

    locs, labels = plt.xticks()
    mylocs = list()
    mylabels = list()

    print(f'{locs}, {labels}')
    for element in coil_set.elements:
        mylocs += (element.position_x - element.length/2, element.position_x, element.position_x + element.length/2)
        mylabels += ('l', element.name, 'r')

    locs = np.concatenate((locs, mylocs))
    labels = np.concatenate((labels, mylabels))

    ax.set_xticks(locs)
    ax.set_xticklabels(labels)
    ax.legend((r'computed $B_x$', r'ideal $B_x$'), loc=1)
    ax.set_xlabel('Neutron Trajectory (m)')
    ax.set_ylabel('Magnetic field (G)')

    # plt.show()
    plt.savefig(f'./bfield_coil_set{plot_name}.png')


def optimize_coils_positions(coil_type, n, max_distance):
    """Compute the optimization.

    Iterate over several distances for each outer coil.


    Parameters
    ----------
    coil_type: Class instance
        Coil Class to be used for the analysis.
    n: int, optional
        Number of iteration steps.
    max_distance: float, optional
        The maximum distance to be iterated over.

    Returns
    -------
    out: np.array
        Array consisting of iteration points for the
    """
    l = define_iteration_values_for_coil_distances(n, max_distance)
    fits = [[0 for i in range(n)] for j in range(n)]

    best_fit = {'fit_value': 0, 'l12': 0, 'l34': 0}

    x_positions = define_computational_grid()

    middle_position = 0

    for i in range(len(l)):
        for j in range(len(l)):
            print(f'compute i= {i} j= {j}')
            # Create CoilSets
            coil_set = CoilSet(coil_type=coil_type,
                               distance_12=l[i], distance_34=l[j],
                               name='CoilSet',
                               position=middle_position)

            # Compute B_field_values
            with Pool(4) as p:
                b_field_values = p.map(coil_set.b_field, x_positions)

            # Ideal b field
            b_max = max(b_field_values)
            ideal_b_field = get_ideal_field(b_max, middle_position, x_positions)

            fit_value = minimizer_function(b_field_values, ideal_b_field)
            fits[j][i] = fit_value

            if fit_value > best_fit['fit_value']:
                best_fit = {'fit_value': fit_value, 'l12': l[i], 'l34': l[j]}

    plot_name = '_test'
    plot_l1l2_cmap(fits, l, plot_name=plot_name)
    plot_ideal_position(middle_position, best_fit, x_positions, coil_type=coil_type, plot_name=plot_name)
