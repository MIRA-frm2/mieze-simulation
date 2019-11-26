# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""The two pairs of coils."""

from elements.base import BasicElement
from elements.coils import Coil

from experiments.mieze.parameters import DISTANCE_2ND_COIL, DISTANCE_3RD_COIL, DISTANCE_4TH_COIL, R_IN, R_OUT, L_IN, L_OUT, N_IN, N_OUT, COIL_SET_CURRENT, absolute_x_position


from utils.helper_functions import save_data_to_file


class CoilSet(BasicElement):
    """Class that implements a coil with more realistic experimental parameters."""

    def __init__(self, position, coil_type=Coil):
        super(CoilSet, self).__init__(position)

        self.coil_type = coil_type

        self._create_coil_set()

    def _create_coil_set(self):
        first_coil_pos = 0

        self.coil_outer_1 = self.coil_type(position=(first_coil_pos, 0, 0), length=L_OUT,
                                           windings=N_OUT, current=-COIL_SET_CURRENT, r=R_OUT, wire_d=0)
        self.coil_inner_1 = self.coil_type(position=(first_coil_pos + DISTANCE_2ND_COIL, 0, 0),
                                           length=L_IN, windings=N_IN, current=COIL_SET_CURRENT, r=R_IN, wire_d=0)
        self.coil_inner_2 = self.coil_type(position=(first_coil_pos + DISTANCE_3RD_COIL, 0, 0),
                                           length=L_IN, windings=N_IN, current=COIL_SET_CURRENT, r=R_IN, wire_d=0)
        self.coil_outer_2 = self.coil_type(position=(first_coil_pos + DISTANCE_4TH_COIL, 0, 0),
                                           length=L_OUT, windings=N_OUT, current=-COIL_SET_CURRENT, r=R_OUT, wire_d=0)

    def b_field(self, x, y, z):
        return self.coil_outer_1.b_field(x, y, z) + self.coil_inner_1.b_field(x, y, z) + self.coil_inner_2.b_field(x, y, z) + self.coil_outer_2.b_field(x, y, z)

    def compute_b_field(self, x_positions):
        b_values = dict()
        for x in absolute_x_position:
            b_values[x, 0, 0] = self.b_field(x, 0, 0)
        return b_values


def compute_default_b_field():
    coil_set = CoilSet(position=0)

    b_field_values = coil_set.compute_b_field(absolute_x_position)
    save_data_to_file(b_field_values, file_name='../data/data')


if __name__ == "__main__":
    compute_default_b_field()
