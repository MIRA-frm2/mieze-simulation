import numpy as np
from scipy.stats import chisquare

from elements.coil_set import CoilSet
from experiments.mieze.parameters import DISTANCE_3RD_COIL, LENGTH_COIL_INNER

import matplotlib.pyplot as plt

from utils.helper_functions import save_data_to_file, unit_square


def minimizer_function(computed_values, expected_values, use_chisquare=False):
    if use_chisquare:
        fit_value = chisquare(1e4 * np.asarray(computed_values),
                              1e4 * np.asarray(expected_values))
    else:
        fit_value = 0
        for j in range(len(computed_values)):
            dif = computed_values[j] - expected_values[j]
            fit_value += dif ** 2
    return fit_value


def optimize_coils_positions():
    n = 2
    max_distance = 0.5

    fits = [[0 for i in range(n)] for j in range(n)]

    l = np.linspace(0, max_distance, n)
    best_fit = {'fit_value': 0, 'l12': 0, 'l34': 0}

    # Computational grid space
    startpoint = -0.5  # [m]
    endpoint = 0.75  # [m]  # Positions.get_position_coilA()
    npoints = 200
    x_positions = np.linspace(startpoint, endpoint, num=npoints)

    for i in range(len(l)):
        for j in range(len(l)):
            print(f'compute i= {i} j= {j}')
            # Create CoilSets
            coil_set = CoilSet(position=0, distance_12=l[i], distance_34=l[j])

            # Compute B_field_values
            b_field_values = coil_set.compute_b_field(x_positions)

            # Ideal b field
            b_max = max(b_field_values)
            ideal_b_field = b_max * unit_square(LENGTH_COIL_INNER/2,
                                                DISTANCE_3RD_COIL - LENGTH_COIL_INNER/2,
                                                x_positions)

            fit_value = minimizer_function(b_field_values, ideal_b_field)
            fits[j][i] = fit_value

            if fit_value > best_fit['fit_value']:
                best_fit = {'fit_value': fit_value, 'l12': l[i], 'l34': l[j]}

    plt.imshow(fits, aspect='auto', extent=[min(l), max(l), min(l), max(l)], origin='lower')
    plt.colorbar()

    plt.show()

    # Plot values
    coil_set = CoilSet(position=0, distance_12=best_fit['l12'], distance_34=best_fit['l34'])

    # Computational grid space
    startpoint = -0.25  # [m]
    endpoint = 0.5  # [m]  # Positions.get_position_coilA()
    npoints = 100
    x_positions = np.linspace(startpoint, endpoint, num=npoints)

    # Compute B_field_values
    b_field_values = coil_set.compute_b_field(x_positions)

    plt.plot(x_positions, b_field_values)

    # Plot idea B field
    b_max = max(b_field_values)
    ideal_b_field = b_max * unit_square(0, DISTANCE_3RD_COIL, x_positions)

    plt.plot(x_positions, ideal_b_field)

    plt.show()


if __name__ == "__main__":
    optimize_coils_positions()
