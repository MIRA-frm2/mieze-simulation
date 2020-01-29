import matplotlib.pyplot as plt
from copy import deepcopy

from simulation.elements.spin_flipper import SpinFlipper
from simulation.experiments.mieze.parameters import absolute_x_position, SPIN_FLIPPER_PARAMETERS


def main():

    kwargs = deepcopy(SPIN_FLIPPER_PARAMETERS)
    kwargs['position'] = 0.5

    spin_flipper = SpinFlipper(**kwargs)

    b_field_x = list()
    b_field_y = list()
    b_field_z = list()

    for val_x in absolute_x_position:
        magnetic_field_vector = spin_flipper.b_field([val_x, 0, 0])

        b_field_x.append(magnetic_field_vector[0])
        b_field_y.append(magnetic_field_vector[1])
        b_field_z.append(magnetic_field_vector[2])

    fig, axs = plt.subplots(ncols=3)

    fig.suptitle('Magnetic fields for the individual element: Spin Flipper')

    axs[0].plot(absolute_x_position, b_field_x)
    axs[0].set_title('Bx')
    axs[1].plot(absolute_x_position, b_field_y)
    axs[1].set_title('By')
    axs[2].plot(absolute_x_position, b_field_z)
    axs[2].set_title('Bz')

    plt.show()


if __name__ == '__main__':
    main()
