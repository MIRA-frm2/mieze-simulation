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

from elements.helmholtz_spin_flipper import HelmholtzSpinFlipper
from elements.coils import Coil

from experiments.mieze.parameters import R_HSF, I_hsf1, HelmholtzSpinFlipper_position_HSF1


class Test(TestCase):

    def setUp(self):
        self.helmholtz_spin_flipper = HelmholtzSpinFlipper(
            coil_type=Coil,
            coil_distance=R_HSF,
            current=I_hsf1,
            position=(HelmholtzSpinFlipper_position_HSF1, 0, 0),
            radius=R_HSF)

    def test_zero_position(self):
        """Test for known reference position values."""

        test_data = {(0, 0, 0): (4.4467801407790776e-05, 0., 0.)}

        for position, reference_value in test_data.items():
            # Evaluate
            b_field_value = self.helmholtz_spin_flipper.b_field(*position)

            # Assert tests
            # Assert data type
            assert isinstance(b_field_value, type(array([0])))
            # Assert numerical value
            self.assertEqual(list(reference_value), list(b_field_value))

    def test_reference_value(self):
        """Test for known reference values."""

        self.helmholtz_spin_flipper.windings = 124
        self.helmholtz_spin_flipper.radius = 14.9 * 1e-2  # [m]
        self.helmholtz_spin_flipper.current = 1

        position = (0, 0, 0)

        reference_value = 7.49 * 1e-4  # [T]

        # Evaluate
        b_field_value = self.helmholtz_spin_flipper.b_field(*position)[0]

        # Assert numerical value
        self.assertEqual(reference_value, b_field_value)

