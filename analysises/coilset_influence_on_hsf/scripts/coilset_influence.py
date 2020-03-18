import matplotlib.pyplot as plt
import numpy as np

from experiments.mieze.main_mieze import compute_magnetic_field_mieze
from experiments.mieze.parameters import HelmholtzSpinFlipper_position_HSF1, WIDTH_CBOX

from utils.helper_functions import read_data_from_file, find_nearest


class MyPlotter:
    """Customized Plotter Class."""

    def __init__(self, data_file='../../data/data_magnetic_field.csv'):
        """

        Parameters
        ----------
        data_file: string
            Location of the data file containing the magnetic field values.
        """
        self.data_file = data_file
        self.x_range, self.y_range, self.z_range, self.bx, self.by, self.bz = read_data_from_file(data_file)

        self.theta_values = None
        self.dtheta_dy = None
        self.b_values = None
        self.adiabatic_condition_values = None
        self.diff_data = None

        self.preadjust_values()

    def preadjust_values(self):
        """Adjust values format for plotting/analysis purposes."""
        self.bx = np.asarray(np.abs(self.bx))

    def compute_influence(self, position=HelmholtzSpinFlipper_position_HSF1):
        index = find_nearest(self.x_range, position)
        # print(self.x_range[index])
        # print(self.bx[index])
        max_val = max(self.bx)
        # print(max_val)
        return self.bx[index]/max_val * 100


def find_edges():
    data_file = "../ideal_position_b_field_values.csv"

    plotter = MyPlotter(data_file=data_file)
    # print(WIDTH_CBOX/2)
    print(plotter.compute_influence(position=WIDTH_CBOX/2))


def find_relative_distance_coilset_sf(iteration_values):

    influence_list = list()

    for val in iteration_values:
        data_file = f'/data/test_data_{int(val*100)}.csv'

        compute_magnetic_field_mieze(filename=data_file, coil_set_distance=val)

        # Using the existing Plotter class, read values and compute values
        plotter = MyPlotter(data_file=data_file)

        influence_list.append(plotter.compute_influence())

    plt.plot(np.asarray(iteration_values) * 100, influence_list)
    plt.xlabel('Distance between right HSF and left outer coil (cm)')
    plt.ylabel(r'Ratio of Coil Set value at HSF over Max Coil Set (%)')

    plt.show()


if __name__ == "__main__":
    iteration_values = np.arange(0.0, 0.2 + 0.025, step=0.025)
    find_relative_distance_coilset_sf(iteration_values)

    # find_edges()
