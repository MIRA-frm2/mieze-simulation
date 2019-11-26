import numpy as np


from elements.coil_set import CoilSet

from plotting_scripts.plot import Plotter

from utils.helper_functions import save_data_to_file


def optimize_coils_positions():
    # Create CoilSets
    coil_set = CoilSet(position=0)

    startpoint = -0.25  # [m]
    endpoint = 0.5  # [m]  # Positions.get_position_coilA()
    npoints = 100
    x_positions = np.linspace(startpoint, endpoint, num=npoints)

    b_field_values = coil_set.compute_b_field(x_positions)
    save_data_to_file(b_field_values, file_name='../data/data')

    plotter = Plotter()

    plotter.plot_field_1d_scalar(component='x', xlabel='Position [m]', ylabel='Magnetic field [a.u.]')


if __name__ == "__main__":
    optimize_coils_positions()
