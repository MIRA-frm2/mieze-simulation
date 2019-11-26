import numpy as np


from elements.coil_set import CoilSet
from experiments.mieze.parameters import DISTANCE_3RD_COIL

from plotting_scripts.plot import Plotter

from utils.helper_functions import save_data_to_file, unit_square


def minimizer_function(computed_values, expected_values):
    fit_value = 0
    for j in range(len(computed_values)):
        dif = computed_values[j] - expected_values[j]
        fit_value += dif ** 2
    return fit_value


def optimize_coils_positions():
    # Create CoilSets
    coil_set = CoilSet(position=0)

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

    fit_value = minimizer_function(b_field_values, ideal_b_field)

    print(fit_value)

    # Plot values
    # save_data_to_file(b_field_values, file_name='../data/data')
    # plotter = Plotter()
    # plotter.plot_field_1d_scalar(component='x', xlabel='Position [m]', ylabel='Magnetic field [a.u.]')


if __name__ == "__main__":
    optimize_coils_positions()
