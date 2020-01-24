# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""A plotting scripts that plots individually the magnetic field for the MIEZE setup elements."""

import matplotlib.pyplot as plt

from simulation.plotting_scripts.plotter import Plotter


def main():
    element_1 = Plotter(filename='../../data/elements_magnetic_fields/data_magnetic_field_Polariser.csv')
    element_2 = Plotter(filename='../../data/elements_magnetic_fields/data_magnetic_field_HelmHoltzSpinFlipper1.csv')
    element_3 = Plotter(filename='../../data/elements_magnetic_fields/data_magnetic_field_SpinFlipper.csv')
    element_4 = Plotter(filename='../../data/elements_magnetic_fields/data_magnetic_field_CoilSet.csv')

    fig, axs = plt.subplots(4, 3)

    fig.suptitle('Magnetic fields for the individual elements')

    axs[0, 0].plot(element_1.x_range, element_1.bx)
    axs[0, 0].set_title('Bx Polariser')
    axs[0, 0].set(ylabel='B [G]')
    axs[0, 1].plot(element_1.x_range, element_1.by)
    axs[0, 1].set_title('By Polariser')
    axs[0, 2].plot(element_1.x_range, element_1.bz)
    axs[0, 2].set_title('Bz Polariser')

    axs[1, 0].plot(element_2.x_range, element_2.bx)
    axs[1, 0].set_title('Bx HSF')
    axs[1, 0].set(ylabel='B [G]')
    axs[1, 1].plot(element_2.x_range, element_2.by)
    axs[1, 1].set_title('By HSF')
    axs[1, 2].plot(element_2.x_range, element_2.bz)
    axs[1, 2].set_title('Bz HSF')

    axs[2, 0].plot(element_3.x_range, element_3.bx)
    axs[2, 0].set_title('Bx SF')
    axs[2, 0].set(ylabel='B [G]')
    axs[2, 1].plot(element_3.x_range, element_3.by)
    axs[2, 1].set_title('By SF')
    axs[2, 2].plot(element_3.x_range, element_3.bz)
    axs[2, 2].set_title('Bz SF')

    axs[3, 0].plot(element_4.x_range, element_4.bx)
    axs[3, 0].set_title('Bx CoilSet')
    axs[3, 0].set(xlabel='x [m]', ylabel='B [G]')
    axs[3, 1].plot(element_4.x_range, element_4.by)
    axs[3, 1].set_title('By CoilSet')
    axs[3, 1].set(xlabel='x [m]', ylabel='B [G]')
    axs[3, 2].plot(element_4.x_range, element_4.bz)
    axs[3, 2].set_title('Bz CoilSet')
    axs[3, 2].set(xlabel='x [m]', ylabel='B [G]')

    plt.show()


if __name__ == "__main__":
    main()
