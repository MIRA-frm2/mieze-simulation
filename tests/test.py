from numpy import sqrt, pi, array

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

    def test_infinity_value(self):
        points = [(0.13, 0.0, 0.0), (-0.13, 0.0, 0.0), (-0.1299999999999999, 0.0, 0.0), (0.13000000000000034, 0.0, 0.0),
                  (0.0885000000000003, 0.0, 0.0), (-0.08849999999999986, 0.0, 0.0)]

        self.square_coil.create_coil(coil_mid_pos=0, length=0.1, windings=1000, current=10,  r=0.05, wire_d=0.006,
                                     angle_y=0, angle_z=0)

        for point in points:
            point = self.square_coil.transform_cylindrical_to_cartesian(*point)
            test_value = self.square_coil.b_field(*point)

            print(test_value)
