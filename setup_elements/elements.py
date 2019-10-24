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
import warnings

from physics_constants import MU_0, pi


class BaseCoil:
    """Class that implements basic method and attributes for a coil."""

    MU_0 = MU_0

    def __init__(self):
        self.prefactor = self.MU_0
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
    def _rotate(vektor, phi, axis):
        """Rotate the vector along an axis for the given angle phi.

        Parameters
        ----------
        vektor: np.array
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
        return np.dot(r, vektor)

    def change_current(self, current):
        """Change the assigned current value."""
        self.current = current

    @staticmethod
    def transform_cartesian_to_cylindrical(x, y, z):
        x = x
        rho = np.sqrt(y ** 2 + z ** 2)
        theta = np.arctan(y / z)
        return x, rho, theta

    @staticmethod
    def transform_cylindrical_to_cartesian(x, rho, theta):
        """Transform coordinates from cylindrical to cartesian."""
        x = x
        y = rho * np.cos(theta)
        z = rho * np.sin(theta)
        return x, y, z

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

        prefactor = 1e4 * self.MU_0 / (2 * np.pi) * self.n * self.current / self.length

        return prefactor * (-self.sum_element_x(n, rho, x, s=-1) + self.sum_element_x(n, rho, x, s=1))

    def sum_element_x(self, n, rho, x, s):
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


class SquareCoil(BaseCoil):

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

        self.prefactor *= 1 / (4 * np.pi) * self.n * self.current / self.length

    def b_field(self, x, rho, theta):
        """Compute the magnetic field."""
        self.x, self.y, self.z = self.transform_cylindrical_to_cartesian(x, rho, theta)

        # self.x -= self.coil_mid_pos

        field = self.prefactor * np.array([self.eta_x, self.eta_y, self.eta_z])

        return self._fix_numerical_error(field)

    @staticmethod
    def _fix_numerical_error(field):
        """Fix computational error.

        Parameters
        ----------
        field: np.array
            3D vector
        """

        # if field[0] > 0:
        #     field[0] = 0
        # elif field[0] < 0:
        #     field[0] = - field[0]

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
    def eta_x(self):
        _eta = 0
        for alpha in np.arange(1, 5):
            alpha_index = alpha - 1

            with warnings.catch_warnings(record=True) as w:
                # Cause all warnings to always be triggered.
                warnings.simplefilter("always")
                # Verify some things
                t1 = ((-1) ** alpha) * self.d[alpha_index] / \
                     (self.r_values[alpha_index] * (self.r_values[alpha_index]
                                                    + (-1) ** (alpha + 1) * self.c[alpha_index]))
                t2 = self.c[alpha_index] \
                     / (self.r_values[alpha_index] * (self.r_values[alpha_index] + self.d[alpha_index]))
                if t1 == 0 or t2 == 0:
                    print(f"0 value at {self.x}, {self.y}, {self.z}")

                _eta += t1 - t2

                if len(w) == 1:
                    # ToDo: what to do when dividing by 0?
                    print(f"Numerical error encountered at {self.x, self.y, self.z}")
                    _eta = 0
        return _eta

    @property
    def eta_y(self):
        _eta = 0
        for alpha in np.arange(1, 5):
            alpha_index = alpha - 1

            _eta += ((-1)**(alpha+1)) * self.z\
                    /(self.r_values[alpha_index] * (self.r_values[alpha_index]
                                                    + ((-1) ** (alpha+1)) * self.c[alpha_index]))
        return _eta

    @property
    def eta_z(self):
        _eta = 0
        for alpha in np.arange(1, 5):
            alpha_index = alpha - 1

            with warnings.catch_warnings(record=True) as w:
                _eta += ((-1) ** (alpha + 1)) * self.z / (
                            self.r_values[alpha_index] * (self.r_values[alpha_index] + self.d[alpha_index]))
                if len(w) == 1:
                    # ToDo: what to do when dividing by 0?
                    print(f"Numerical error encountered at {self.x, self.y, self.z}")
                    _eta = 0

        return _eta

    @property
    def c(self):
        return np.array([self.a + self.x,
                         self.a - self.x,
                         -self.a + self.x,
                         -self.a - self.x])

    @property
    def d(self):
        return np.array([self.b + self.y,
                         self.b + self.y,
                         self.y - self.b,
                         self.y - self.b])

    @property
    def r_values(self):
        return np.array([
            np.sqrt((self.a + self.x) ** 2 + (self.b + self.y) ** 2 + self.z ** 2),
            np.sqrt((self.a - self.x) ** 2 + (self.b + self.y) ** 2 + self.z ** 2),
            np.sqrt((self.a - self.x) ** 2 + (self.b - self.y) ** 2 + self.z ** 2),
            np.sqrt((self.a + self.x) ** 2 + (self.b - self.y) ** 2 + self.z ** 2)]
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

    def b_field_rho(self, x, r, rho):
        """ # radial magnetic field

        """
        x -= self.coil_mid_pos
        if rho == 0:
            return 0

        rho = abs(rho)

        return self.prefactor * (self.sum_element_rho(rho, x, s=1) - self.sum_element_rho(rho, x, s=-1))

    def sum_element_rho(self, rho, x, s):
        return self.beta(rho, x, s) / rho \
               * ((2 - self.m(rho, x, s)) * self._k(self.m(rho, x, s)) - 2 * self._e(self.m(rho, x, s)))

    def b_field_x(self, x, r, rho=0):
        """
        # N = turns of winding, I = current (A), R = radius (m), l = length (m)
        # x für mich y

        """
        x -= self.coil_mid_pos

        rho = abs(rho)  # symmetric

        n = 4.0 * r * rho / (rho + r) ** 2

        prefac = 1e4 * self.MU_0 / (2 * np.pi) * self.n * self.current / self.length

        return prefac * (-self.sum_element_x(n, rho, x, s=-1) + self.sum_element_x(n, rho, x, s=1))

    def sum_element_x(self, n, rho, x, s):
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
