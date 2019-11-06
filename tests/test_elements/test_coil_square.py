# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Numerical tests for the codebase."""

from numpy import sqrt, pi
from unittest import TestCase

from elements.coils import SquareCoil

from utils.physics_constants import MU_0


class TestReferenceValues(TestCase):

    def setUp(self) -> None:
        self.square_coil = SquareCoil(coil_mid_pos=0, length=1, windings=1, wire_d=0.006)

    def test_square_coil_center_value(self):
        """Test for a reference value in the center of a square coil."""

        # Setup
        current = 10
        side_length = 1

        reference_value = sqrt(2) * MU_0 * current / (pi * side_length)

        numerical_error_acceptance = 1e-8

        self.square_coil.current = current
        self.square_coil.r = side_length/2

        # Evaluate
        test_value = self.square_coil.b_field(0, 0, 0)

        # print(test_value, reference_value)

        # Assert
        assert abs(reference_value - test_value[0]) < numerical_error_acceptance
        assert test_value[1] == 0
        assert test_value[2] == 0
