# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Individual magnetic elements (coils) for the setup."""

import json
import logging
import mpmath
import numpy as np
import warnings

from simulation.elements.base import BasicElement

from utils.helper_functions import get_phi, adjust_field, sanitize_output
from utils.physics_constants import MU_0, pi, factor_T_to_G

# Set the warning filter to errors such that one can catch them as they were errors
# This is used to handle divisions by numerical 0
warnings.filterwarnings('error')

# Create a custom logger
logger = logging.getLogger(__name__)


class BaseCoil(BasicElement):
    """Class that implements basic method and attributes for a coil."""

    def __init__(self, position, name, **kwargs):
        """

        Parameters
        ----------
        Keyword Arguments:
            current: float
                The value of the coil current.
            length: float
                Coil length
            windings: int
                Number of coil windings.
            wire_d: float
                Diameter of the coil wire.
            r_eff: float
                The effective coil radius.
            r_min: float
                The inner (min) radius value for the coil.
            r_max: float
                The outer (max) radius value for the coil.
        """
        super(BaseCoil, self).__init__(position, name)

        self.iteration = 0
        self.prefactor = MU_0

        self.length = kwargs.get('length', None)
        self.current = kwargs.get('current', 0)

        # Place the simplified coil at the middle of the real coil
        # if self.length:
        #     self.position_x += self.length / 2

        self.windings = kwargs.get('windings', None)
        self.wire_d = kwargs.get('wire_d', None)
        self.wire_spacing = kwargs.get('wire_spacing', None)
        self.radial_layers = kwargs.get('radial_layers', None)

        # Compute the effective radius, from possibly different inputs
        r_eff = kwargs.get('r_eff', None)
        r_min = kwargs.get('r_min', None)
        r_max = kwargs.get('r_max', None)

        if r_eff:
            self.r = r_eff
        # # Compute the effective radius from the min, max radiuses and the wire diameter.
        # elif r_min and r_max:
        #     inverse_r_sum = 0
        #     if self.wire_d:
        #         radius_values = np.arange(r_min, r_max, self.wire_d)
        #         num_layers = len(radius_values)
        #         for r_i in radius_values:
        #             inverse_r_sum += 1 / r_i
        #
        #         self.r = 1 / (inverse_r_sum / num_layers)

        elif r_min and r_max:
            self.compute_effective_radius(r_min, r_max)

        else:
            self.width = kwargs.get('width', None)
            self.height = kwargs.get('height', None)

            if not self.width or not self.height:
                raise Exception('Radius value not set.')

        self.angle_y = kwargs.get('angle_y', 0)
        self.angle_z = kwargs.get('angle_z', 0)

    def __repr__(self):
        return json.dumps(self, default=lambda o: o.__dict__,
                          sort_keys=True, indent=4)

    def compute_effective_radius(self, r_min, r_max):
        """Compute the effective radius from the other physical parameters."""
        inverse_r_sum = 0

        num_layers = (r_max - r_min) / (self.wire_d + self.wire_spacing)
        step = (r_max - r_min) / num_layers
        # print(f'r_min:{r_min}\nwire_d:{self.wire_d}\nwire_spacing:{self.wire_spacing}')
        for r_i in np.arange(r_min, r_max, step):
            inverse_r_sum += 1 / r_i

        self.r = 1 / (inverse_r_sum / num_layers)
        # print(self.r)

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
        """Compute zeta function."""
        return x - s * self.length / 2.0

    def beta(self, rho, x, s):
        """Compute beta function."""
        return np.sqrt((rho + self.r) ** 2 + self.zeta(x=x, s=s) ** 2)

    def m(self, rho, x, s):
        """Compute m function."""
        return 4.0 * self.r * rho / self.beta(rho, x, s) ** 2

    def b_field(self, x, y, z):
        raise NotImplementedError


class Coil(BaseCoil):
    """Class that implements an ideal circular coil."""

    name = 'Coil'

    def __init__(self, name, position, **kwargs):
        """Simulate physical geometry of a simplified coil."""

        super(Coil, self).__init__(position, name, **kwargs)

        self.prefactor *= self.windings * self.current / (2 * pi * self.length)

    @sanitize_output
    def b_field_rho(self, x, rho):
        """Computes radial magnetic field

        Parameters
        ----------
        x: float, int
            Position on the x axis
        rho: float, int
            Position in the radial direction (assuming symmetry)

        Returns
        -------
        out: float
            Magnetic field value in radial direction.
        """
        if rho == 0:
            return 0

        x -= self.position_x
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

    @sanitize_output
    def b_field_x(self, x, rho=0):
        """Compute the magnetic field in x direction.

        Parameters
        ----------
        x: int, float
            Position along the cylinder.
        rho: int, float
            Position in the radial direction.

        Returns
        -------

        """
        x -= self.position_x

        # symmetric
        rho = abs(rho)

        # print(f'r:{self.r}\nrho:{rho}')
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
        """Compute the magnetic field given the position in cartesian coordinates."""
        r = np.sqrt(z ** 2 + y ** 2)

        phi = get_phi(y, z)
        # print(f'y:{z}\ny:{y}\nr:{r}')
        field = np.array((self.b_field_x(x, r),
                          self.b_field_rho(x, r) * np.cos(phi),
                          self.b_field_rho(x, r) * np.sin(phi))
                         )

        if self.angle_y != 0:
            field = self._rotate(field, self.angle_y, np.array([0, 1, 0]))

        if self.angle_z != 0:
            field = self._rotate(field, self.angle_z, np.array([0, 0, 1]))

        return factor_T_to_G * field


class RealCoil(Coil):
    """Class that implements a coil with more realistic experimental parameters."""

    def b_field_rho(self, x, rho):
        """Compute the magnetic field in rho (radial) direction.


        Parameters
        ----------
        x: float
            Value on the x axis, along the coil.
        rho: float
            Radius value.
        """
        if rho == 0:
            return 0

        x -= self.position_x

        return self.prefactor * (self.sum_element_rho(rho, x, s=1) - self.sum_element_rho(rho, x, s=-1))

    def sum_element_rho(self, rho, x, s):
        """Compute the summation for the rho direction."""
        return self.beta(rho, x, s) \
            / rho * ((2 - self.m(rho, x, s)) * self._k(self.m(rho, x, s)) - 2 * self._e(self.m(rho, x, s)))

    def b_field_beam(self, x, r, rho=0):
        """Compute the magnetic field in the x direction.

        # N = turns of winding, I = current (A), R = radius (m), l = length (m)
        # x für mich y

        """
        x -= self.position_x

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
            * ((rho - self.r) / (rho + self.r) * self._p(n, self.m(rho, x, s)) - self._k(self.m(rho, x, s)))\
               

    def b_field(self, x, y, z):
        """Compute the magnetic field given the position in cartesian coordinates."""
        self.iteration += 1

        field = np.array([0., 0., 0.])
        r = np.sqrt(y ** 2 + z ** 2)
        phi = get_phi(y, z)

        x_positions = np.arange(-(self.windings * self.wire_d / 2 - self.wire_d / 2),
                                self.windings * self.wire_d / 2,
                                self.wire_d)

        rs = np.arange(self.r, self.r + self.wire_d * (self.windings - 0.5), self.wire_d)

        for x0 in x_positions:
            for R in rs:
                print(f'iteration:{self.iteration} with x:{x0} and r:{R}')
                # logger.error(f'{field}\n{field_add}')
                field_add = np.array([self.b_field_beam(x + x0, R, r),
                                      self.b_field_rho(x + x0, R) * np.cos(phi),
                                      self.b_field_rho(x + x0, R) * np.sin(phi)])
                # logger.error(f'{field}\n{field_add}')
                field += field_add
                # field = np.add(field, field_add, out=field, casting='unsafe')

        if self.angle_y != 0:
            field = self._rotate(field, self.angle_y, np.array([0, 1, 0]))

        if self.angle_z != 0:
            field = self._rotate(field, self.angle_z, np.array([0, 0, 1]))

        return factor_T_to_G * field


class RectangularCoil(BaseCoil):
    """Class that implements a rectangular coil."""

    def __init__(self, position, name, **kwargs):
        """Simulate physical geometry of the coil."""

        super(RectangularCoil, self).__init__(position, name, **kwargs)

        self.x, self.y, self.z = None, None, None

        self.prefactor *= 1 / (4 * np.pi) * self.windings * self.current / self.length

    def check_physical_coil_overlap(self):
        """Deal with division by 0.

        If the magnetic field has to be numerically computed at the position of the coil, this implies a division by
        0 leading to numerical errors. This function is used to handle such cases by setting the magnetic field to 0.
        """
        # print(self.x, self.y, self.z, self.width, self.height)
        numerical_error_acceptance = min(self.width, self.height) * 1e-2

        if abs(self.y - self.width) < numerical_error_acceptance \
                or abs(self.y + self.width) < numerical_error_acceptance:
            return True

        if abs(self.x - self.height) < numerical_error_acceptance\
                or abs(self.x + self.height) < numerical_error_acceptance:
            return True

    @staticmethod
    def _reverse_coordinates(x, y, z):
        """Fix coordinates for the equations that have the coordinate axes reversed."""
        return np.array([z, y, x])

    def b_field(self, x, y, z):
        """Compute the magnetic field.

        Returns
        -------
        magnetic field in cartesian coordinates
        """
        # self.x, self.y, self.z = x, y, z
        self.x, self.y, self.z = self._reverse_coordinates(x, y, z)

        field = self.check_physical_coil_overlap()
        if field:
            return np.array([0, 0, 0])

        self.x -= self.position_x

        field = self.prefactor * np.array([self.sum_element_x, self.sum_element_y, self.sum_element_z])

        field = adjust_field(field)

        return factor_T_to_G * self._reverse_coordinates(*self.sanitize_output(field))

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
            return np.array([0, 0, 0])
        if any(field) > 1e6 or any(field) < -1e6:
            return np.array([0, 0, 0])

        return field

    @property
    def sum_element_z(self):
        """Compute the summation for the z direction.

        Returns
        -------
        out: float
        """
        _eta = 0
        for alpha in np.arange(1, 5):
            # Python indexing starts from 0, so whenever alpha is used as an array index it has to be subtracted by 1
            alpha_index = alpha - 1

            try:
                t1 = ((-1) ** alpha) * self.d[alpha_index] / \
                     (self.r_values[alpha_index] * (self.r_values[alpha_index]
                                                    + (-1) ** (alpha + 1) * self.c[alpha_index]))
            except RuntimeWarning:
                t1 = 0

            try:
                t2 = self.c[alpha_index] \
                     / (self.r_values[alpha_index] * (self.r_values[alpha_index] + self.d[alpha_index]))
            except RuntimeWarning:
                t2 = 0

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

            try:
                _eta += ((-1)**(alpha+1)) * self.z / (self.r_values[alpha_index] * (self.r_values[alpha_index]
                                                      + ((-1) ** (alpha+1)) * self.c[alpha_index]))
            except RuntimeWarning:
                _eta = 0

        return _eta

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

            try:
                _eta += ((-1) ** (alpha + 1)) * self.z / (
                            self.r_values[alpha_index] * (self.r_values[alpha_index] + self.d[alpha_index]))
            except RuntimeWarning:
                _eta = 0
        return _eta

    @property
    def c(self):
        """Compute height parameter from the actual coil"""
        return np.array([self.height + self.x,
                         self.height - self.x,
                         -self.height + self.x,
                         -self.height - self.x])

    @property
    def d(self):
        """Compute width parameter from the actual coil."""
        return np.array([self.y + self.width,
                         self.y + self.width,
                         self.y - self.width,
                         self.y - self.width])

    @property
    def r_values(self):
        """Compute the distance from the actual coil."""
        return np.array([
            np.sqrt((self.width + self.y) ** 2 + (self.height + self.x) ** 2 + self.z ** 2),
            np.sqrt((self.width + self.y) ** 2 + (self.height - self.x) ** 2 + self.z ** 2),
            np.sqrt((self.width - self.y) ** 2 + (self.height - self.x) ** 2 + self.z ** 2),
            np.sqrt((self.width - self.y) ** 2 + (self.height + self.x) ** 2 + self.z ** 2)]
        )
