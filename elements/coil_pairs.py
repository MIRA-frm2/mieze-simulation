# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""The two pairs of coils."""

import numpy as np


from elements.base import BasicElement
from elements.coils import Coil

from experiments.mieze.parameters import DISTANCE_2ND_COIL, DISTANCE_3RD_COIL, DISTANCE_4TH_COIL, R_IN, R_OUT, L_IN, L_OUT, N_IN, N_OUT, COIL_SET_CURRENT


from plotting_scripts.plot import Plotter

from utils.helper_functions import save_data_to_file


class CoilSet(BasicElement):
    """Class that implements a coil with more realistic experimental parameters."""

    def __init__(self, position, coil_type=Coil):
        super(CoilSet, self).__init__(position)

        self.coil_type = coil_type

        self.first_coil_pos = self.position_x

        self.elements = list()

        self._create_coil_set()

    def _create_coil_set(self):
        self._create_coil_inner_set()
        self._create_coil_outer_set()

    def _create_coil_inner_set(self):
        coil_inner_1 = self.coil_type(position=(self.first_coil_pos + DISTANCE_2ND_COIL, 0, 0),
                                           length=L_IN, windings=N_IN, current=COIL_SET_CURRENT, r=R_IN, wire_d=0)
        coil_inner_2 = self.coil_type(position=(self.first_coil_pos + DISTANCE_3RD_COIL, 0, 0),
                                           length=L_IN, windings=N_IN, current=COIL_SET_CURRENT, r=R_IN, wire_d=0)

        self.elements.append(coil_inner_1)
        self.elements.append(coil_inner_2)

    def _create_coil_outer_set(self):
        coil_outer_1 = self.coil_type(position=(self.first_coil_pos, 0, 0), length=L_OUT,
                                           windings=N_OUT, current=-COIL_SET_CURRENT, r=R_OUT, wire_d=0)
        coil_outer_2 = self.coil_type(position=(self.first_coil_pos + DISTANCE_4TH_COIL, 0, 0),
                                           length=L_OUT, windings=N_OUT, current=-COIL_SET_CURRENT, r=R_OUT, wire_d=0)

        self.elements.append(coil_outer_1)
        self.elements.append(coil_outer_2)

    def b_field(self, x, y, z):
        b_field = 0
        for element in self.elements:
            b_field += element.b_field(x, y, z)
        return b_field

    def compute_b_field(self, x_positions):
        b_values = dict()
        for x in x_positions:
            b_values[x, 0, 0] = self.b_field(x, 0, 0)
        return b_values


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
