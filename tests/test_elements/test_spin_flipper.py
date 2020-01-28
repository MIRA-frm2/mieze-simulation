# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Numerical tests for the codebase."""

from unittest import TestCase

from simulation.elements.spin_flipper import SpinFlipper


class Test(TestCase):

    def setUp(self) -> None:
        self.spin_flipper = SpinFlipper(name='TestSpinFlipper', current=1, r_eff=1, windings=1, length=1, width=1, height=1)

    def test_zero_position(self):
        """Test for known reference values."""

        test_data = {1.0: -0.00023420086393421903}

        for position, reference_value in test_data.items():
            # Evaluate
            b_field_value = self.spin_flipper.b_field([position, 0, 0])
            # print(f'{b_field_value} {reference_value}')

            # Assert numerical value
            self.assertEqual(reference_value, b_field_value[0])
