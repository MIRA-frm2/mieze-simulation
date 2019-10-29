# -*- coding: utf-8 -*-
#
# This file is part of the E21 FRM2 Research Group.
# Copyright (C) 2019 TUM.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the License; see LICENSE file for more details.


"""Individual magnetic elements (coils) for the setup."""

import mpmath
import numpy as np

from setup_elements.physics_constants import MU_0
from setup_elements.helper_functions import transform_cylindrical_to_cartesian


class BaseCoil:
    """Class that implements basic method and attributes for a coil."""

    def __init__(self):
        self.prefactor = MU_0
        self.coil_mid_pos = None
        self.length = None
        self.n = None
        self.windings_x = None
        self.windings_y = None

        self.r = None
        self.current = None
        self.wire_d = None
        self.angle_y = None
        self.angle_z = None

        self.a = None
        self.b = None

        self.x = None
        self.y = None
        self.z = None

    def create_coil(self):
        """Simulate physical geometry of the coil."""
        raise NotImplementedError

    @staticmethod
    def _k(m):
        """Compute the elliptical function of first kind."""
        return float(mpmath.ellipk(m))

    @staticmethod
    def _e(m):
        """Compute the elliptical function of second kind."""
        return float(mpmath.ellipe(m))

    @staticmethod
    def _p(n, m):
        """Compute the elliptical function third kind"""
        return float(mpmath.ellippi(n, np.pi / 2, m))

    @staticmethod
    def _rotate(vector, phi, axis):
        """Rotate the vector along an axis for the given angle phi.

        Parameters
        ----------
        vector: np.array
            3D vector
        phi: float, int
            Angle in radians
        axis: np.array
            3D vector used as reference for the rotation

        Returns
        -------
        out: np.array
            The rotated vector
        """
        n = axis / np.linalg.norm(axis)
        c = np.cos(phi)
        s = np.sin(phi)
        n1 = n[0]
        n2 = n[1]
        n3 = n[2]
        r = [[n1 ** 2 * (1 - c) + c, n1 * n2 * (1 - c) - n3 * s, n1 * n3 * (1 - c) + n2 * s],
             [n2 * n1 * (1 - c) + n3 * s, n2 ** 2 * (1 - c) + c, n2 * n3 * (1 - c) - n1 * s],
             [n3 * n1 * (1 - c) - n2 * s, n3 * n2 * (1 - c) + n1 * s, n3 ** 2 * (1 - c) + c]]
        return np.dot(r, vector)

    def change_current(self, current):
        """Change the assigned current value."""
        self.current = current

    def zeta(self, x, s):
        return x - s * self.length / 2.0

    def beta(self, rho, x, s):
        return np.sqrt((rho + self.r) ** 2 + self.zeta(x=x, s=s) ** 2)

    def m(self, rho, x, s):
        return 4.0 * self.r * rho / self.beta(rho, x, s) ** 2


class Coil(BaseCoil):
    """Class that implements an ideal circular coil."""

    def create_coil(self, coil_mid_pos=0, length=0.1, windings=100, current=10, r=0.05, wire_d=0.006, angle_y=0,
                    angle_z=0):
        """Simulate physical geometry of the coil."""
        self.coil_mid_pos = coil_mid_pos
        self.length = length
        self.n = windings

        inverse_r_sum = 0
        num_layers = windings * wire_d / length
        r_large = r + wire_d * num_layers

        for r_i in np.arange(r + wire_d / 2, r_large, wire_d):
            inverse_r_sum += 1 / r_i
        self.r = 1 / (inverse_r_sum / num_layers)
        self.current = current
        self.angle_y = angle_y
        self.angle_z = angle_z

        self.prefactor *= 1e4 / (2 * np.pi) * self.n * self.current / self.length

    def b_field_rho(self, x, rho):
        """Computes radial magnetic field

        Parameters
        ----------
        x: float, int
            Position on the x axis
        rho: float, int
            Position

        Returns
        -------

        """
        x -= self.coil_mid_pos
        if rho == 0:
            return 0

        rho = abs(rho)

        return self.prefactor * (self.sum_element_rho(rho, x, s=1) - self.sum_element_rho(rho, x, s=-1))

    def sum_element_rho(self, rho, x, s):
        """Compute the summation for the rho direction.

        Returns
        -------
        out: float
        """
        return self.beta(rho, x, s) / rho \
            * ((2 - self.m(rho, x, s)) * self._k(self.m(rho, x, s)) - 2 * self._e(self.m(rho, x, s)))

    def b_field_x(self, x, rho=0):
        """Compute the magnetic field in x direction.


        # N = turns of winding, I = current (A), R = radius (m), l = length (m)

        Parameters
        ----------
        x: int, float

        rho: int, float
            Position along the cylinder.

        Returns
        -------

        """
        x -= self.coil_mid_pos

        rho = abs(rho)  # symmetric

        n = 4.0 * self.r * rho / (rho + self.r) ** 2

        return self.prefactor * (-self.sum_element_x(n, rho, x, s=-1) + self.sum_element_x(n, rho, x, s=1))

    def sum_element_x(self, n, rho, x, s):
        """Compute the summation for the x direction.

        Returns
        -------
        out: float
        """
        return self.zeta(x, s) / self.beta(rho, x, s) \
            * ((rho - self.r) / (rho + self.r) * self._p(n, self.m(rho, x, s)) - self._k(self.m(rho, x, s)))

    def b_field(self, x, y, z):
        """Compute the magnetic field given the position."""
        r = np.sqrt(y ** 2 + z ** 2)
        field = np.array((self.b_field_x(x, r), self.b_field_rho(x, y), self.b_field_rho(x, z)))

        if self.angle_y != 0:
            field = self._rotate(field, self.angle_y, np.array([0, 1, 0]))

        if self.angle_z != 0:
            field = self._rotate(field, self.angle_z, np.array([0, 0, 1]))

        return field
        # return transform_cylindrical_to_cartesian(*field)


class SquareCoil(BaseCoil):
    """Class that implements a square coil.

    The coordinates from the reference paper are changed from (x, y, z) to (y, z, x), hence the change in the formulae.
    """

    def create_coil(self, coil_mid_pos=0, length=0.1, windings=100, current=10, r=0.05, wire_d=0.006, angle_y=0,
                    angle_z=0):
        """Simulate physical geometry of the coil."""
        self.coil_mid_pos = coil_mid_pos

        self.current = current

        self.length = length

        self.a = r
        self.b = r

        self.n = windings

        self.r = r

        self.angle_y = angle_y
        self.angle_z = angle_z

        self.prefactor *= 1e4 / (4 * np.pi) * self.n * self.current / self.length

    def b_field(self, x, rho, theta):
        """Compute the magnetic field.

        Parameters
        ----------
        position in cylindrical coordinates

        Returns
        -------
        magnetic field in cartesian coordinates
        """
        self.x, self.y, self.z = transform_cylindrical_to_cartesian(x, rho, theta)

        self.x -= self.coil_mid_pos

        field = self.prefactor * np.array([self.sum_element_x, self.sum_element_y, self.sum_element_z])

        return self.sanitize_output(field)

    @staticmethod
    def sanitize_output(field):
        """Sanitize magnetic field..

        Remove possible numerical singularities.

        Parameters
        ----------
        field: np.array
            3D vector
        """

        if any(np.isinf(field)) or any(np.isnan(field)):
            field = np.array([0, 0, 0])

        return field

    # @property
    # def b_field_x(self):
    #     """Compute the B magnetic field in x direction."""
    #     return self.eta_x
    #
    # @property
    # def b_field_y(self):
    #     """Compute the B magnetic field in y direction."""
    #     return self.eta_y
    #
    # @property
    # def b_field_z(self):
    #     """Compute the B magnetic field in z direction."""
    #     return  self.eta_z

    @property
    def sum_element_x(self):
        """Compute the summation for the x direction.

        Returns
        -------
        out: float
        """
        _eta = 0
        for alpha in np.arange(1, 5):
            alpha_index = alpha - 1

            t1 = ((-1) ** alpha) * self.d[alpha_index] / \
                 (self.r_values[alpha_index] * (self.r_values[alpha_index]
                                                + (-1) ** (alpha + 1) * self.c[alpha_index]))
            t2 = self.c[alpha_index] / (self.r_values[alpha_index] * (self.r_values[alpha_index] + self.d[alpha_index]))

            _eta += t1 - t2

        return _eta

    @property
    def sum_element_y(self):
        """Compute the summation for the y direction.

        Returns
        -------
        out: float
        """
        _eta = 0
        for alpha in np.arange(1, 5):
            alpha_index = alpha - 1

            _eta += ((-1)**(alpha+1)) * self.z / (self.r_values[alpha_index] * (self.r_values[alpha_index]
                                                  + ((-1) ** (alpha+1)) * self.c[alpha_index]))
        return _eta

    @property
    def sum_element_z(self):
        """Compute the summation for the z direction.

        Returns
        -------
        out: float
        """
        _eta = 0
        for alpha in np.arange(1, 5):
            alpha_index = alpha - 1

            _eta += ((-1) ** (alpha + 1)) * self.z / (
                        self.r_values[alpha_index] * (self.r_values[alpha_index] + self.d[alpha_index]))

        return _eta

    @property
    def c(self):
        return np.array([self.a + self.y,
                         self.a - self.y,
                         -self.a + self.y,
                         -self.a - self.y])

    @property
    def d(self):
        return np.array([self.b + self.z,
                         self.b + self.z,
                         self.z - self.b,
                         self.z - self.b])

    @property
    def r_values(self):
        return np.array([
            np.sqrt((self.a + self.y) ** 2 + (self.b + self.z) ** 2 + self.x ** 2),
            np.sqrt((self.a - self.y) ** 2 + (self.b + self.z) ** 2 + self.x ** 2),
            np.sqrt((self.a - self.y) ** 2 + (self.b - self.z) ** 2 + self.x ** 2),
            np.sqrt((self.a + self.y) ** 2 + (self.b - self.z) ** 2 + self.x ** 2)]
        )


class RealCoil(BaseCoil):
    """Class that implements a coil with more realistic experimental parameters."""
    MU_0 = 4e-7 * np.pi

    # oc = Oct2Py()
    def create_coil(self, coil_mid_pos=0, length=0.1, windings=100, current=10, r=0.05, wire_d=0.006, angle_y=0,
                    angle_z=0):
        """Simulate physical geometry of the coil."""
        self.coil_mid_pos = coil_mid_pos
        self.length = length
        self.n = windings
        self.windings_x = round(length / wire_d)
        self.windings_y = windings / self.windings_x
        self.r = r
        self.current = current
        self.wire_d = wire_d
        self.angle_y = angle_y
        self.angle_z = angle_z

        self.prefactor *= 1e4 / (2 * np.pi) * self.n * self.current / self.length

    def b_field_rho(self, x, r, axis_point):
        """Compute the magnetic field in rho (radial) direction.


        Parameters
        ----------
        x: float
            Value on the x axis, along the coil.
        r: float
            Radius value.
        axis_point: float
            Value on the y or z axis.
        """
        x -= self.coil_mid_pos
        if axis_point == 0:
            return 0

        if r:
            # ToDo: is the radius really not playing a role?
            pass

        rho = abs(axis_point)

        return self.prefactor * (self.sum_element_rho(rho, x, s=1) - self.sum_element_rho(rho, x, s=-1))

    def sum_element_rho(self, rho, x, s):
        """Compute the summation for the rho direction."""
        return self.beta(rho, x, s) \
            / rho * ((2 - self.m(rho, x, s)) * self._k(self.m(rho, x, s)) - 2 * self._e(self.m(rho, x, s)))

    def b_field_x(self, x, r, rho=0):
        """Compute the magnetic field in the x direction.

        # N = turns of winding, I = current (A), R = radius (m), l = length (m)
        # x für mich y

        """
        x -= self.coil_mid_pos

        rho = abs(rho)  # symmetric

        n = 4.0 * r * rho / (rho + r) ** 2

        return self.prefactor * (-self.sum_element_x(n, rho, x, s=-1) + self.sum_element_x(n, rho, x, s=1))

    def sum_element_x(self, n, rho, x, s):
        """Compute the summation for the x direction.

        Parameters
        ----------
        n: int
            Winding number
        rho: float
            Radius value.
        x: float
            Position on the x axis along the coil.
        s: [-1, 1]
            Boundary values for hte elliptical functions
        """
        return self.zeta(x, s) / self.beta(rho, x, s) \
            * ((rho - self.r) / (rho + self.r) * self._p(n, self.m(rho, x, s)) - self._k(self.m(rho, x, s)))

    def b_field(self, x, y, z):
        """Compute the magnetic field given the position."""
        field = np.array((0., 0., 0.))
        r = np.sqrt(y ** 2 + z ** 2)

        # self.pbar.update(1)
        x_positions = np.arange(-(self.windings_x * self.wire_d / 2 - self.wire_d / 2),
                                self.windings_x * self.wire_d / 2, self.wire_d)
        # print(len(x_positions))
        rs = np.arange(self.r, self.r + self.wire_d * (self.windings_y - 0.5), self.wire_d)
        # print(len(rs))
        for x0 in x_positions:
            for R in rs:
                field_add = np.array((self.b_field_x(x + x0, R, r),
                                      self.b_field_rho(x + x0, R, y),
                                      self.b_field_rho(x + x0, R, z)))
                field = np.add(field, field_add, out=field)

        if self.angle_y != 0:
            field = self._rotate(field, self.angle_y, np.array([0, 1, 0]))

        if self.angle_z != 0:
            field = self._rotate(field, self.angle_z, np.array([0, 0, 1]))

        return field
