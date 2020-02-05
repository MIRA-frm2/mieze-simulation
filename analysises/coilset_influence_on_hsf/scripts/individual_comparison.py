import matplotlib.pyplot as plt
import numpy as np

from experiments.mieze.parameters import ELEMENTS_POSITIONS_ABSOLUTE

from utils.helper_functions import read_data_from_file, find_nearest


def plot(x, y):
    plt.plot(x * 100, y * 100, '*-')

    plt.xlabel('Distance between centre HSF and centre CoilSet [cm]')
    plt.ylabel('Ratio of the Tail Magnetic Field from CoilSet to Max Value of HSF [%]')

    plt.show()


def main():
    folder_location = '../../../data/elements_magnetic_fields'

    distance_between_centers = ELEMENTS_POSITIONS_ABSOLUTE['coil_set'] - ELEMENTS_POSITIONS_ABSOLUTE['spin_flipper']

    x, y, z, bx, by, bz = read_data_from_file(f'{folder_location}/data_magnetic_field_HelmHoltzSpinFlipper1.csv')
    index_sf = find_nearest(x, ELEMENTS_POSITIONS_ABSOLUTE['spin_flipper'])
    b_spin_flipper = bx[index_sf]

    x, y, z, bx_coilset, by, bz = read_data_from_file(f'{folder_location}/data_magnetic_field_CoilSet.csv')

    print(distance_between_centers)
    d = x[:index_sf]
    d.reverse()
    d = distance_between_centers + np.asarray(d)
    y = np.asarray(bx[:index_sf]) / b_spin_flipper

    plot(d, y)


if __name__ == "__main__":
    main()
