import numpy as np


from elements.coil_set import CoilSet
from experiments.mieze.parameters import DISTANCE_3RD_COIL

import matplotlib.pyplot as plt

from utils.helper_functions import save_data_to_file, unit_square


def minimizer_function(computed_values, expected_values):
    fit_value = 0
    for j in range(len(computed_values)):
        dif = computed_values[j] - expected_values[j]
        fit_value += dif ** 2
    return fit_value


def optimize_coils_positions():
    n = 25
    max_distance = 0.25

    fits = [[0 for i in range(n)] for j in range(n)]

    l = np.linspace(0, max_distance, n)

    for i in range(len(l)):
        for j in range(len(l)):
            print(f'compute i= {i} j= {j}')
            # Create CoilSets
            coil_set = CoilSet(position=0, distance_2nd_coil=l[i], distance_4th_coil=l[j])

            # Computational grid space
            startpoint = -0.25  # [m]
            endpoint = 0.5  # [m]  # Positions.get_position_coilA()
            npoints = 100
            x_positions = np.linspace(startpoint, endpoint, num=npoints)

            # Compute B_field_values
            b_field_values = coil_set.compute_b_field(x_positions)

            # Ideal b field
            b_max = max(b_field_values)
            ideal_b_field = unit_square(0, DISTANCE_3RD_COIL, x_positions)

            fit_value = minimizer_function(b_field_values, ideal_b_field * b_max)
            fits[j][i] = fit_value

    plt.imshow(fits, aspect='auto')  #, extent=[min(x), max(x), min(y), max(y)])
    plt.colorbar()

    plt.show()

    # Plot values
    # save_data_to_file(b_field_values, file_name='../data/data')
    # plotter = Plotter()
    # plotter.plot_field_1d_scalar(component='x', xlabel='Position [m]', ylabel='Magnetic field [a.u.]')


if __name__ == "__main__":
    optimize_coils_positions()
