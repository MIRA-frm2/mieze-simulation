# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Investigation script to find the ideal position for the two Coil Set Pairs."""

import matplotlib.pyplot as plt
from multiprocessing import Pool
import numpy as np
from scipy.stats import chisquare

from elements.coil_set import CoilSet
from experiments.mieze.parameters import DISTANCE_2ND_COIL, DISTANCE_3RD_COIL, LENGTH_COIL_INNER, DISTANCE_BETWEEN_INNER_COILS


from utils.helper_functions import save_data_to_file, unit_square


def minimizer_function(computed_values, expected_values, use_chisquare=False):
    """Evaluate the similarity of the computed and expected arrays.

    Parameters
    ----------
    computed_values: np.array
        The array of computed values.
    expected_values: np.array
        The array of ideal/expeted values, used as reference.
    use_chisquare: bool
        Flag indicating whether to use the chisquare test.
        In order to use chisquare, one needs arrays with at least a value of 5 in each array element.

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
        for j in range(len(computed_values)):
            dif = computed_values[j] - expected_values[j]
            fit_value += dif ** 2
        fit_value = 1 / fit_value

    return fit_value


def define_iteration_values(n=1):
    if n == 1:
        l = (0.0,)
    else:
        max_distance = 0.1
        step = max_distance / n

        l = np.linspace(0, max_distance, n)
        # l = np.linspace(0, 1, 1)

    return l


def define_computational_grid():
    # Computational grid space
    startpoint = -0.25  # [m]
    endpoint = 0.25  # [m]  # Positions.get_position_coilA()
    return np.linspace(startpoint, endpoint, num=200)


def plot_l1l2_cmap(fits, l):
    if len(l) > 1:
        plt.imshow(fits, aspect='auto', extent=[min(l), max(l), min(l), max(l)], origin='lower')
        plt.colorbar()
        plt.xlabel('L1 (First Outer Coil Distance) [m]')
        plt.ylabel('L2 (Second Outer Coil Distance) [m]')

        # plt.show()
        plt.savefig('../docs/experiments/MIEZE/coil_set_optimization_fine.png')
        plt.close()


def plot_ideal_position(middle_position, best_fit, x_positions):
    # Plot values
    print(best_fit)
    coil_set = CoilSet(name='CoilSet', position=middle_position, distance_12=best_fit['l12'], distance_34=best_fit['l34'])

    # Compute B_field_values
    with Pool(4) as p:
        b_field_values = p.map(coil_set.b_field, x_positions)

    print(b_field_values)

    fig = plt.figure()
    ax = plt.axes()

    plt.plot(x_positions, b_field_values)

    # Plot idea B field
    b_max = max(b_field_values)
    ideal_b_field = get_ideal_field(b_max, middle_position, x_positions)

    plt.plot(x_positions, ideal_b_field)

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
    plt.savefig('../docs/experiments/MIEZE/bfield_coil_set_fine.png')


def get_ideal_field(b_max, middle_position, x_positions):
    return b_max * unit_square(middle_position - DISTANCE_BETWEEN_INNER_COILS / 2 - LENGTH_COIL_INNER - LENGTH_COIL_INNER * 0.5,
                               middle_position + DISTANCE_BETWEEN_INNER_COILS / 2 + LENGTH_COIL_INNER * 0.5,
                               x_positions)


def optimize_coils_positions():
    """Compute the optimization.

    Iterate over several distances for each outer coil.
    """
    n = 10
    l = define_iteration_values(n)
    fits = [[0 for i in range(n)] for j in range(n)]

    best_fit = {'fit_value': 0, 'l12': 0, 'l34': 0}

    x_positions = define_computational_grid()

    middle_position = 0

    for i in range(len(l)):
        for j in range(len(l)):
            print(f'compute i= {i} j= {j}')
            # Create CoilSets
            coil_set = CoilSet(distance_12=l[i], distance_34=l[j], name='CoilSet', position=middle_position)

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

    plot_l1l2_cmap(fits, l)
    plot_ideal_position(middle_position, best_fit, x_positions)


if __name__ == "__main__":
    optimize_coils_positions()
