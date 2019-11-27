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

from experiments.mieze.parameters import DISTANCE_2ND_COIL, DISTANCE_3RD_COIL, DISTANCE_4TH_COIL, R_IN, R_OUT, L_IN, L_OUT, N_IN, N_OUT, COIL_SET_CURRENT


class CoilSet(BasicElement):
    """Class that implements a coil with more realistic experimental parameters."""

    def __init__(self, position, coil_type=Coil, distance_12=None, distance_34=None):
        super(CoilSet, self).__init__(position)

        self.coil_type = coil_type

        self.first_inner_coil_pos = self.position_x

        self.distance_12 = - distance_12 if distance_12 else (- DISTANCE_2ND_COIL)
        self.distance_23 = DISTANCE_3RD_COIL - DISTANCE_2ND_COIL

        self.distance_34 = distance_34 if distance_34 else (DISTANCE_4TH_COIL - DISTANCE_3RD_COIL)
        self.distance_34 += self.distance_23

        self.elements = list()

        self._create_coil_set()

    def _create_coil_set(self):
        self._create_coil_inner_set()
        self._create_coil_outer_set()

    def _create_coil_inner_set(self):
        coil_inner_1 = self.coil_type(position=self.first_inner_coil_pos,
                                      length=L_IN, windings=N_IN, current=COIL_SET_CURRENT, r=R_IN, wire_d=0)
        coil_inner_2 = self.coil_type(position=self.first_inner_coil_pos + (DISTANCE_3RD_COIL-DISTANCE_2ND_COIL),
                                      length=L_IN, windings=N_IN, current=COIL_SET_CURRENT, r=R_IN, wire_d=0)

        self.elements.append(coil_inner_1)
        self.elements.append(coil_inner_2)

    def _create_coil_outer_set(self):
        coil_outer_1 = self.coil_type(position=self.first_inner_coil_pos + self.distance_12, length=L_OUT,
                                      windings=N_OUT, current=-COIL_SET_CURRENT, r=R_OUT, wire_d=0)
        coil_outer_2 = self.coil_type(position=self.first_inner_coil_pos + self.distance_34,
                                      length=L_OUT, windings=N_OUT, current=-COIL_SET_CURRENT, r=R_OUT, wire_d=0)

        self.elements.append(coil_outer_1)
        self.elements.append(coil_outer_2)

    def b_field(self, x, y, z):
        b_field = 0
        for element in self.elements:
            b_field += element.b_field(x, y, z)[0]
        return b_field

    def compute_b_field(self, x_positions):
        b_values = list()
        for x in x_positions:
            b_values.append(self.b_field(x, 0, 0))
        return b_values
