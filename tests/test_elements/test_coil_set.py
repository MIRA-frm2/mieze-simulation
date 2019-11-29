import matplotlib.pyplot as plt
import numpy as np

from elements.coil_set import CoilSet


def simple_execution():
    coil_set = CoilSet(name='CoilSet', position=0, distance_12=0.5, distance_34=0.5)

    # Computational grid space
    startpoint = -0.75  # [m]
    endpoint = 0.75  # [m]  # Positions.get_position_coilA()
    npoints = 100
    x_positions = np.linspace(startpoint, endpoint, num=npoints)

    # Compute B_field_values
    b_field_values = coil_set.compute_b_field(x_positions)

    plt.plot(x_positions, b_field_values)
    plt.show()


if __name__ == "__main__":
    simple_execution()
