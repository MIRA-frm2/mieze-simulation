# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""The two pairs of coils."""

from elements.base import BasicElement
from elements.coils import Coil, RealCoil

from experiments.mieze.parameters import DISTANCE_2ND_COIL, DISTANCE_3RD_COIL, DISTANCE_4TH_COIL, RADIUS_COIL_INNER_EFFECTIVE, RADIUS_COIL_OUTER_EFFECTIVE, LENGTH_COIL_INNER, LENGTH_COIL_OUTER, N_WINDINGS_COIL_INNER, N_WINDINGS_COIL_OUTER, COIL_SET_CURRENT, RADIUS_COIL_INNER_MIN, RADIUS_COIL_INNER_MAX, WIRE_D


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
        """Create the two coil pairs."""
        self._create_coil_inner_set()
        self._create_coil_outer_set()

    def _create_coil_inner_set(self):
        """Create the two inner coil pairs."""
        coil_inner_1 = self.coil_type(current=COIL_SET_CURRENT,
                                      length=LENGTH_COIL_INNER,
                                      position=self.first_inner_coil_pos,
                                      r_min=RADIUS_COIL_INNER_MIN, r_max=RADIUS_COIL_INNER_MAX,
                                      windings=N_WINDINGS_COIL_INNER,
                                      wire_d=WIRE_D)
        coil_inner_2 = self.coil_type(current=COIL_SET_CURRENT,
                                      length=LENGTH_COIL_INNER,
                                      position=self.first_inner_coil_pos + (DISTANCE_3RD_COIL-DISTANCE_2ND_COIL),
                                      r_min=RADIUS_COIL_INNER_MIN, r_max=RADIUS_COIL_INNER_MAX,
                                      windings=N_WINDINGS_COIL_INNER,
                                      wire_d=WIRE_D)

        self.elements.append(coil_inner_1)
        self.elements.append(coil_inner_2)

    def _create_coil_outer_set(self):
        """Create the two outer coil pairs."""
        coil_outer_1 = self.coil_type(current=-COIL_SET_CURRENT,
                                      length=LENGTH_COIL_OUTER,
                                      position=-LENGTH_COIL_OUTER + self.first_inner_coil_pos + self.distance_12,
                                      r_eff=RADIUS_COIL_OUTER_EFFECTIVE,
                                      windings=N_WINDINGS_COIL_OUTER,
                                      wire_d=WIRE_D)

        coil_outer_2 = self.coil_type(current=-COIL_SET_CURRENT,
                                      length=LENGTH_COIL_OUTER,
                                      position=LENGTH_COIL_OUTER+self.first_inner_coil_pos + self.distance_34,
                                      r_eff=RADIUS_COIL_OUTER_EFFECTIVE,
                                      windings=N_WINDINGS_COIL_OUTER,
                                      wire_d=WIRE_D)

        self.elements.append(coil_outer_1)
        self.elements.append(coil_outer_2)

    def b_field(self, x, y, z):
        """Compute the magnetic field for one specific position.

        Parameters
        ----------
        x
        y
        z

        Returns
        -------

        """
        b_field = 0
        for element in self.elements:
            b_field += element.b_field(x, y, z)[0]
        return b_field

    def compute_b_field(self, x_positions):
        """Compute the magnetic field for the array of positions."""
        b_values = list()
        for x in x_positions:
            b_values.append(self.b_field(x, 0, 0))
        return b_values
