from numpy import sqrt, pi

from unittest import TestCase

from setup_elements.elements import SquareCoil
from physics_constants import MU_0


class TestReferenceValues(TestCase):

    def setUp(self) -> None:
        self.square_coil = SquareCoil()

    def test_square_coil_center_value(self):
        """Test for a reference value in the center of a square coil."""
        current = 10
        side_length = 1

        reference_value = sqrt(2) * MU_0 * current / (pi * side_length)

        numerical_error_acceptance = 1e-5

        self.square_coil.create_coil(coil_mid_pos=0, length=1, windings=1, current=current, r=side_length, wire_d=0.006)
        test_value = self.square_coil.b_field(0, 0, 0)

        # print(test_value, reference_value)

        assert abs(reference_value - test_value[0]) < numerical_error_acceptance
        assert test_value[1] == 0
        assert test_value[2] == 0
