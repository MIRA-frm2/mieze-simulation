# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Numerical tests for the codebase."""

from numpy import array
from unittest import TestCase

from elements.polariser import Polariser


class Test(TestCase):

    def setUp(self) -> None:
        self.polariser = Polariser()

    def test_zero_position(self):
        """Test for known reference values."""

        test_data = {(0, 0, 0): (0., -719.9999999999998, 0.)}

        for position, reference_value in test_data.items():
            # Evaluate
            b_field_value = self.polariser.b_field(*position)
            # print(b_field_value)

            # Assert tests
            # Assert data type
            assert isinstance(b_field_value, type(array([0])))
            # Assert numerical value
            self.assertEqual(list(reference_value), list(b_field_value))

