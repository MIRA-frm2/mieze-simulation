# -*- coding: utf-8 -*-
#
# This file is part of the E21 FRM2 Research Group.
# Copyright (C) 2019 TUM.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.


"""Description """

import mpmath
import numpy as np

from physics_constants import MU_0, pi


class BaseCoil:
    """Class that implements basic method and attributes for a coil."""

    def __init__(self):
        self.coil_mid_pos = None
        self.l = None
        self.N = None
        self.windings_x = None
        self.windings_y = None
        self.R = None
        self.I = None
        self.wire_d = None
        self.angle_y = None
        self.angle_z = None

    def create_coil(self):
        """Simulate physical geometry of the coil."""
        raise NotImplementedError

    def _K(cls, m):
        """Compute the elliptical function of first kind."""
        return float(mpmath.ellipk(m))

    def _E(cls, m):
        """Compute the elliptical function of second kind."""
        return float(mpmath.ellipe(m))

    def _P(cls, n, m):
        """Compute the elliptical function third kind"""
        return float(mpmath.ellippi(n, np.pi / 2, m))

    def _rotate(self, vektor, phi, axis):
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
        R = [[n1 ** 2 * (1 - c) + c, n1 * n2 * (1 - c) - n3 * s, n1 * n3 * (1 - c) + n2 * s],
             [n2 * n1 * (1 - c) + n3 * s, n2 ** 2 * (1 - c) + c, n2 * n3 * (1 - c) - n1 * s],
             [n3 * n1 * (1 - c) - n2 * s, n3 * n2 * (1 - c) + n1 * s, n3 ** 2 * (1 - c) + c]]
        return np.dot(R, vektor)

    def change_current(self, I):
        """Change the assigned current value."""
        self.I = I


class Coil(BaseCoil):
    """Class that implements an ideal circular coil."""

    MU_0 = 4e-7 * np.pi

    # oc = Oct2Py()

    def create_coil(self, coil_mid_pos=0, length=0.1, windings=100, current=10, r=0.05, wire_d=0.006, angle_y=0,
                    angle_z=0):
        """Simulate physical geometry of the coil."""
        self.coil_mid_pos = coil_mid_pos
        self.l = length
        self.N = windings
        inverse_r_sum = 0
        num_layers = windings * wire_d / length
        R = r + wire_d * num_layers

        for r_i in np.arange(r + wire_d / 2, R, wire_d):
            inverse_r_sum += 1 / r_i

        self.R = 1 / (inverse_r_sum / num_layers)
        self.I = current
        self.angle_y = angle_y
        self.angle_z = angle_z


    def Brho(self, x, rho):
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
        beta = lambda s: np.sqrt((rho + self.R) ** 2 + (x - s * self.l / 2) ** 2)
        m = lambda s: 4 * self.R * rho / beta(s) ** 2
        sum_element = lambda s: beta(s) / rho * ((2 - m(s)) * self._K(m(s)) - 2 * self._E(m(s)))
        prefac = 1e4 * self.MU_0 / (2 * np.pi) * self.N * self.I / self.l

        return prefac * (sum_element(1) - sum_element(-1))

    def Bx(self, x, rho=0):
        """Compute the magnetic field in x direction.

        # N = turns of winding, I = current (A), R = radius (m), l = length (m)
        # x für mich y

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
        zeta = lambda s: x - s * self.l / 2.0

        beta = lambda s: np.sqrt((rho + self.R) ** 2 + zeta(s) ** 2)

        m = lambda s: 4.0 * self.R * rho / beta(s) ** 2

        n = 4.0 * self.R * rho / (rho + self.R) ** 2

        prefac = 1e4 * self.MU_0 / (2 * np.pi) * self.N * self.I / self.l

        sum_element = lambda s: zeta(s) / beta(s) * ((rho - self.R) / (rho + self.R) * self._P(n, m(s)) - self._K(m(s)))

        return prefac * (-sum_element(-1) + sum_element(1))

    def B_field(self, x, y, z):
        """Compute the magnetic field given the position."""
        r = np.sqrt(y ** 2 + z ** 2)
        field = np.array((self.Bx(x, r), self.Brho(x, y), self.Brho(x, z)))
        if self.angle_y != 0:
            field = self._rotate(field, self.angle_y, np.array([0, 1, 0]))

        if self.angle_z != 0:
            field = self._rotate(field, self.angle_z, np.array([0, 0, 1]))

        return field


class SquareCoil:

    prefactor = MU_0 / (4 * pi)

    def __init__(self):
        self.a = None
        self.b = None

        self.I = None

    def create_coil(self, coil_mid_pos, length, windings, current,  r, wire_d, angle_y, angle_z):
        """Simulate physical geometry of the coil."""
        self.I = current

        self.a = length
        self.b = length

    def B_field(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        return self.Bx, self.By, self.By

    @property
    def Bx(self):
        return self.prefactor * self.eta_x

    @property
    def By(self):
        return 0

    @property
    def Bz(self):
        return self.prefactor * self.I * self.eta_z

    @property
    def eta_x(self):
        _eta = 0
        for alpha in np.arange(1, 5):
            _eta += ((-1)**(alpha+1)) * self.z/(self.r[alpha-1] * (self.r[alpha-1] + self.d[alpha-1]))
        return _eta

    @property
    def eta_z(self):
        _eta = 0
        for alpha in np.arange(1, 5):
            t1 = ((-1) ** alpha) * self.d[alpha-1] / (self.r[alpha-1] * (self.r[alpha-1] + (-1)**(alpha +1) * self.C[alpha-1]))
            t2 = self.C[alpha-1] / (self.r[alpha-1] * (self.r[alpha-1] + self.d[alpha-1]))
            _eta += t1 - t2
        return _eta

    @property
    def C(self):
        return np.array([self.a+self.x, self.a-self.x, -self.a + self.x, -self.a - self.x])

    @property
    def d(self):
        return np.array([self.b + self.y, self.b + self.y, self.y - self.b, self.y - self.b])

    @property
    def r(self):
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
        self.l = length
        self.windings_x = round(length / wire_d)
        self.windings_y = windings / self.windings_x
        self.R = r
        self.I = current
        self.wire_d = wire_d
        self.angle_y = angle_y
        self.angle_z = angle_z

    def Brho(self, x, R, rho):  # radial magnetic field
        x -= self.coil_mid_pos
        if rho == 0:
            return 0

        rho = abs(rho)
        beta = lambda s: np.sqrt((rho + R) ** 2 + (x - s * self.l / 2) ** 2)
        m = lambda s: 4 * R * rho / beta(s) ** 2
        sum_element = lambda s: beta(s) / rho * ((2 - m(s)) * self._K(m(s)) - 2 * self._E(m(s)))
        prefac = 1e4 * self.MU_0 / (2 * np.pi) * self.N * self.I / self.l

        return prefac * (sum_element(1) - sum_element(-1))

    def Bx(self, x, R, rho=0):  # N = turns of winding, I = current (A), R = radius (m), l = length (m)
        # x für mich y
        x -= self.coil_mid_pos

        rho = abs(rho)  # symmetric
        zeta = lambda s: x - s * self.l / 2.0

        beta = lambda s: np.sqrt((rho + R) ** 2 + zeta(s) ** 2)

        m = lambda s: 4.0 * R * rho / beta(s) ** 2

        n = 4.0 * R * rho / (rho + R) ** 2

        prefac = 1e4 * self.MU_0 / (2 * np.pi) * self.N * self.I / self.l

        sum_element = lambda s: zeta(s) / beta(s) * ((rho - R) / (rho + R) * self._P(n, m(s)) - self._K(m(s)))

        return prefac * (-sum_element(-1) + sum_element(1))

    def B_field(self, x, y, z):
        """Compute the magnetic field given the position."""
        field = np.array((0., 0., 0.))
        r = np.sqrt(y ** 2 + z ** 2)

        # self.pbar.update(1)
        x_positions = np.arange(-(self.windings_x * self.wire_d / 2 - self.wire_d / 2),
                                self.windings_x * self.wire_d / 2, self.wire_d)
        print(len(x_positions))
        Rs = np.arange(self.R, self.R + self.wire_d * (self.windings_y - 0.5), self.wire_d)
        print(len(Rs))
        for x0 in x_positions:
            for R in Rs:
                field_add = np.array((self.Bx(x + x0, R, r), self.Brho(x + x0, R, y), self.Brho(x + x0, R, z)))
                field = np.add(field, field_add, out=field)

        if self.angle_y != 0:
            field = self._rotate(field, self.angle_y, np.array([0, 1, 0]))

        if self.angle_z != 0:
            field = self._rotate(field, self.angle_z, np.array([0, 0, 1]))

        return field

