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

from experiments.mieze.parameters import LENGTH_COIL_INNER, LENGTH_COIL_OUTER, N_WINDINGS_COIL_INNER, N_WINDINGS_COIL_OUTER, COIL_SET_CURRENT, RADIUS_COIL_INNER_MIN, RADIUS_COIL_INNER_MAX, WIRE_D, DISTANCE_BETWEEN_INNER_COILS, RADIAL_LAYERS, WIRE_SPACING, RADIUS_COIL_OUTER_MIN, RADIUS_COIL_OUTER_MAX


class CoilSet(BasicElement):
    """Class that implements a coil with more realistic experimental parameters."""

    def __init__(self, name, position, coil_type=Coil, distance_12=None, distance_34=None):
        super(CoilSet, self).__init__(position, name)

        self.coil_type = coil_type

        self.middle = self.position_x

        min_distance_inner_outer_coils = (LENGTH_COIL_INNER + LENGTH_COIL_OUTER)/2
        self.distance_12 = min_distance_inner_outer_coils + distance_12

        self.distance_34 = min_distance_inner_outer_coils + distance_34

        self.elements = list()

        self._create_coil_set()

    def _create_coil_set(self):
        """Create the two coil pairs."""
        self._create_coil_inner_set()
        self._create_coil_outer_set()

    def compute_coil_position(self, name):
        """Compute the relative position of each coil."""
        if name == 'coil_inner_1':
            return self.middle - DISTANCE_BETWEEN_INNER_COILS / 2 - LENGTH_COIL_INNER
        elif name == 'coil_inner_2':
            return self.middle + DISTANCE_BETWEEN_INNER_COILS / 2
        elif name == 'coil_outer_1':
            try:
                return self.coil_inner_1.position_x - self.distance_12
            except AttributeError:
                return -LENGTH_COIL_OUTER + self.distance_12
        elif name == 'coil_outer_2':
            try:
                return self.coil_inner_2.position_x + self.distance_34
            except AttributeError:
                return LENGTH_COIL_OUTER + self.distance_34

    def _create_coil_inner_set(self):
        """Create the two inner coil pairs."""
        self.coil_inner_1 = self.coil_type(current=COIL_SET_CURRENT,
                                           length=LENGTH_COIL_INNER,
                                           name='I1',
                                           position=self.compute_coil_position('coil_inner_1'),
                                           r_min=RADIUS_COIL_INNER_MIN, r_max=RADIUS_COIL_INNER_MAX,
                                           radial_layers=RADIAL_LAYERS,
                                           windings=N_WINDINGS_COIL_INNER,
                                           wire_d=WIRE_D,
                                           wire_spacing=WIRE_SPACING)
        self.coil_inner_2 = self.coil_type(current=COIL_SET_CURRENT,
                                           length=LENGTH_COIL_INNER,
                                           name='I2',
                                           position=self.compute_coil_position('coil_inner_2'),
                                           radial_layers=RADIAL_LAYERS,
                                           r_min=RADIUS_COIL_INNER_MIN,
                                           r_max=RADIUS_COIL_INNER_MAX,
                                           windings=N_WINDINGS_COIL_INNER,
                                           wire_d=WIRE_D,
                                           wire_spacing=WIRE_SPACING)

        self.elements.append(self.coil_inner_1)
        self.elements.append(self.coil_inner_2)

    def _create_coil_outer_set(self):
        """Create the two outer coil pairs."""

        coil_outer_1 = self.coil_type(current=-COIL_SET_CURRENT,
                                      length=LENGTH_COIL_OUTER,
                                      name='O1',
                                      position=self.compute_coil_position('coil_outer_1'),
                                      r_min=RADIUS_COIL_OUTER_MIN,
                                      r_max=RADIUS_COIL_OUTER_MAX,
                                      windings=N_WINDINGS_COIL_OUTER,
                                      wire_d=WIRE_D,
                                      wire_spacing=WIRE_SPACING)

        coil_outer_2 = self.coil_type(current=-COIL_SET_CURRENT,
                                      length=LENGTH_COIL_OUTER,
                                      name='O2',
                                      position=self.compute_coil_position('coil_outer_2'),
                                      r_min=RADIUS_COIL_OUTER_MIN,
                                      r_max=RADIUS_COIL_OUTER_MAX,
                                      windings=N_WINDINGS_COIL_OUTER,
                                      wire_d=WIRE_D,
                                      wire_spacing=WIRE_SPACING)

        self.elements.append(coil_outer_1)
        self.elements.append(coil_outer_2)

    def b_field(self, x, y=0, z=0):
        """Compute the magnetic field for one specific position.

        Parameters
        ----------
        x
        y
        z

        Returns
        -------
        b_field: float
            The value of the magnetic field along x direction.
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
