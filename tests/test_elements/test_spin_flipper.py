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

from elements.spin_flipper import SpinFlipper


class Test(TestCase):

    def setUp(self) -> None:
        self.spin_flipper = SpinFlipper()
        self.spin_flipper.create_element()

    def test_zero_position(self):
        """Test for known reference values."""

        test_data = {1.0: -0.0019916348542965423}

        for position, reference_value in test_data.items():
            # Evaluate
            b_field_value = self.spin_flipper.sf(sf_name='sf1', z_position=position)
            # print(b_field_value)

            # Assert numerical value
            self.assertEqual(reference_value, b_field_value)

