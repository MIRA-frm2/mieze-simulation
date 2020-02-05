# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019, 2020 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""A plotting scripts that plots the magnetic field one element."""

import matplotlib.pyplot as plt
from copy import deepcopy

from simulation.elements.helmholtz_pair import HelmholtzPair
from simulation.elements.spin_flipper import SpinFlipper
from experiments.mieze.parameters import absolute_x_position, SPIN_FLIPPER_PARAMETERS, HELMHOLTZCOILS_PARAMETERS


def compute_element_magnetic_field(element_class, init_kwargs):
    """Compute the magnetic field."""
    spin_flipper = element_class(**init_kwargs)

    b_field_x = list()
    b_field_y = list()
    b_field_z = list()

    # Compute the magnetic field at each position.
    for val_x in absolute_x_position:
        magnetic_field_vector = spin_flipper.b_field([val_x, 0, 0])

        b_field_x.append(magnetic_field_vector[0])
        b_field_y.append(magnetic_field_vector[1])
        b_field_z.append(magnetic_field_vector[2])

    return b_field_x, b_field_y, b_field_z


def plot_magnetic_field(element_name, b_field_x, b_field_y, b_field_z):
    """Plot the given magnetic field."""

    fig, axs = plt.subplots(ncols=3)

    # Plot title
    fig.suptitle(f'Magnetic fields for the individual element: {element_name}')

    axs[0].plot(absolute_x_position, b_field_x)
    axs[0].set_title('$B_x$')
    axs[1].plot(absolute_x_position, b_field_y)
    axs[1].set_title('$B_y$')
    axs[2].plot(absolute_x_position, b_field_z)
    axs[2].set_title('$B_z$')

    # Set x and y axis for the plots
    axs[0].set(ylabel='B [G]')

    for ax in axs:
        ax.set(xlabel='x [m]')

    plt.show()


def plot_element_magnetic_field(element_class, init_kwargs):
    """Plot the magnetic field for the given element class.

    Parameters
    ----------
    element_class: class
    init_kwargs: dict
        Dictionary containing keyword arguments for the class constructor.
    """
    b_field_x, b_field_y, b_field_z = compute_element_magnetic_field(element_class, init_kwargs)
    plot_magnetic_field(element_class.class_name, b_field_x, b_field_y, b_field_z)


def plot_spin_flipper():
    """Plot the magnetic field for the Spin Flipper."""
    init_kwargs = deepcopy(SPIN_FLIPPER_PARAMETERS)
    init_kwargs['position'] = 0.5

    plot_element_magnetic_field(element_class=SpinFlipper, init_kwargs=init_kwargs)


def plot_helmholtz_coils():
    """Plot the magnetic field for the Helmholtz Coils."""
    init_kwargs = deepcopy(HELMHOLTZCOILS_PARAMETERS)
    # init_kwargs['position'] = 0.5

    plot_element_magnetic_field(element_class=HelmholtzPair, init_kwargs=init_kwargs)


if __name__ == '__main__':
    # plot_spin_flipper()
    plot_helmholtz_coils()
