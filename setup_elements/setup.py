# -*- coding: utf-8 -*-
#
# This file is part of the E21 FRM2 Research Group.
# Copyright (C) 2019 TUM.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the License; see LICENSE file for more details.


"""General setup allowing placement of different elements."""

import itertools
import matplotlib.pyplot as plt
from multiprocessing import Pool
import numpy as np
from scipy.interpolate import interp1d

from setup_elements.elements import Coil, RealCoil, SquareCoil


class Setup:
    """Class that simulates a physical setup."""
    def __init__(self, increment, coil_type='simple'):
        self.elements = []
        self.b = None
        self.setup_changed = False
        self.x_range = None
        self.y_range = None
        self.z_range = None
        self.current = None
        self.increment = increment
        self.coil_type = coil_type

    def create_coil(self, coil_mid_pos=0, length=0.1, windings=1000, current=10,  r=0.05, wire_d=0.006, angle_y=0,
                    angle_z=0):
        """Create the physical geometry of the coils."""
        if self.coil_type == 'simple':
            element = Coil()
        elif self.coil_type == 'square':
            element = SquareCoil()
        elif self.coil_type == 'real':
            element = RealCoil()
        else:
            raise Exception('Coil type not recognized.')

        element.create_coil(coil_mid_pos, length, windings, current,  r, wire_d, angle_y, angle_z)
        self.elements.append(element)
        self.setup_changed = True

    @staticmethod
    def _find_nearest(array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return array[idx]

    def b_x(self, x, rho=0):
        """Compute magnetic field in x direction."""
        field = 0
        for element in self.elements:
            field += element.b_x(x, rho)
        return field

    def b_rho(self, x, rho):
        """Compute magnetic field in rho direction."""
        field = 0
        for element in self.elements:
            np.add(field, element.b_rho(x, rho), out=field)

        return field

    def b_field(self, r: '(x, y, z)'):
        """Compute magnetic field."""
        field = np.array((0., 0., 0.))
        for element in self.elements:
            np.add(field, element.b_field(*r), out=field)
            # field += element.B_field(*r)

        # print(f'The total computed magnetic field at point {r} is: {field}')
        return field

    def change_current(self, current):
        """Change the current valaue."""
        self.current = current
        self.setup_changed = True

    def calculate_b_field(self, zero, meshsize: '(x, y, z)'):
        """Calculate the magnetic field."""
        
        if self.b is None or self.setup_changed:

            self.x_range = np.arange(zero, meshsize[0] + zero, self.increment)
            self.y_range = np.arange(-meshsize[1], meshsize[1]+self.increment, self.increment)
            self.z_range = np.arange(-meshsize[2], meshsize[2]+self.increment, self.increment)
            # print(self.x_range, self.y_range, self.z_range)
            args = list(itertools.product(self.x_range, self.y_range, self.z_range))
            print(len(args), " calculations")
            print('calculate')

            with Pool(4) as p:
                result = p.map(self.b_field, args)

            print('calculation finished')

            self.b = dict(zip(args, result))

        if meshsize[0]-1 > max(self.x_range) or meshsize[1]-1 > max(self.y_range) or meshsize[2]-1 > max(self.z_range):

            self.x_range = np.arange(zero, meshsize[0] + zero, self.increment)
            self.y_range = np.arange(-meshsize[1], meshsize[1], self.increment)
            self.z_range = np.arange(-meshsize[2], meshsize[2], self.increment)
            args = list(itertools.product(self.x_range, self.y_range, self.z_range))

            for key in self.b.keys():
                if key in args:
                    args.remove(key)

            print('calculate extension')

            with Pool(4) as p:
                result = p.map(self.b_field, args)

            print('calculation finished')

            self.b.update(dict(zip(args, result)))

        self.setup_changed = False
    #
    # def return_1d_fi(self, start, end, rho=0):
    #     self.calculate_b_field(zero=start, meshsize=(end - start, rho, 0))
    #
    #     x_index_start = np.where(self.x_range == self._find_nearest(self.x_range, start))
    #     x_index_end = np.where(self.x_range == self._find_nearest(self.x_range, end))
    #
    #     x_values = self.x_range[int(x_index_start[0]):int(x_index_end[0])+1]
    #
    #     if rho != 0:
    #         y_values = [np.linalg.norm(self.b[(x, rho, 0)]) for x in x_values]
    #     else:
    #         y_values = [self.b[(x, 0, 0)][0] for x in x_values]
    #
    #     return sum(y_values)*self.increment

    def plot_1d_abs(self, start, end, rho=0):

        self.calculate_b_field(zero=start, meshsize=(end - start, rho, 0))

        x_index_start = np.where(self.x_range == self._find_nearest(self.x_range, start))
        x_index_end = np.where(self.x_range == self._find_nearest(self.x_range, end))

        x_values = self.x_range[int(x_index_start[0]):int(x_index_end[0]) + 1]

        if rho != 0:
            y_values = [np.linalg.norm(self.get_b((x, rho, 0))) for x in x_values]
        else:
            try:
                y_values = [self.b[(x, 0, 0)][0] for x in x_values]
            except KeyError:
                # Todo: How to plot b field for non-zero y and z values?
                y_values = None
        plt.plot(self.x_range, y_values)
        plt.show()

    def plot_1d_vector(self, start, end, rho=0):

        self.calculate_b_field(zero=start, meshsize=(end - start, rho, 0))
        x_index_start = np.where(self.x_range == self._find_nearest(self.x_range, start))
        x_index_end = np.where(self.x_range == self._find_nearest(self.x_range, end))

        x_pos = self.x_range[int(x_index_start[0]):int(x_index_end[0]) + 1]

        y_pos = 0
        b_vec = np.array([self.b[(x, rho, 0)][:2] for x in x_pos])
        u = b_vec[:, 0]
        v = b_vec[:, 1]
        print(u, v)

        u_plot_data = np.array([0] * len(u))
        v_plot_data = np.array([0] * len(v))
        for i in range(len(u)):
            magnitude = np.sqrt(u[i] ** 2 + v[i] ** 2)
            if magnitude == 0:
                u_plot_data[i] = 0
                v_plot_data[i] = 0
            else:
                u_plot_data[i] = u[i]/magnitude
                v_plot_data[i] = v[i]/magnitude
        plt.quiver(x_pos, y_pos, u_plot_data, v_plot_data)
        plt.show()

    def data_compare_1d(self, x_data, y_data, rho=0):

        self.calculate_b_field(zero=min(x_data), meshsize=(max(x_data), rho, 0))

        x_values = self.x_range
        # print(x_data, x_values)

        if rho != 0:
            y_values = [np.linalg.norm(self.b[(x, rho, 0)]) for x in x_values]

        else:
            y_values = [self.b[(x, rho, 0)][0] for x in x_values]

        b = interp1d(x_data, y_data)

        b_data = np.array([b(x) for x in x_values])
        plt.plot(x_values, y_values-b_data)
        plt.show()

    # def plot_2d_vector(self, start, end, rho=0):
    #
    #     self.calculate_b_field(zero=start, meshsize=(end - start, rho, 0))
    #
    #     x_index_start = np.where(self.x_range == self._find_nearest(self.x_range, start))
    #     x_index_end = np.where(self.x_range == self._find_nearest(self.x_range, end))
    #
    #     x_value = self.x_range[int(x_index_start[0]):int(x_index_end[0]) + 1]
    #     distance_of_arrows = int(len(self.y_range)/25.)+1
    #
    #     y_pos = self.y_range[::distance_of_arrows]
    #
    #     x_pos = x_value[0::int(1.*len(x_value)/len(y_pos))+1]
    #
    #     b_vec = np.array([[self.b[(x, y, 0)] for x in x_pos] for y in y_pos])
    #     print(x_pos)
    #     print(y_pos)
    #     print(self.b, b_vec)
    #
    #     u = b_vec[:, :, 0]
    #     v = b_vec[:, :, 1]
    #
    #     plt.quiver(x_pos, y_pos, u/np.sqrt(u**2+v**2), v/np.sqrt(u**2+v**2))
    #     plt.show()
    #
    # def plot_2d_map(self, start, end, rho):
    #
    #     self.calculate_b_field(zero=start, meshsize=(end - start, rho, 0))
    #     x_index_start = np.where(self.x_range == self._find_nearest(self.x_range, start))
    #     x_index_end = np.where(self.x_range == self._find_nearest(self.x_range, end))
    #
    #     x_value = self.x_range[int(x_index_start[0]):int(x_index_end[0]) + 1]
    #
    #     # y_index_start = np.where(self.y_range == self.find_nearest(self.y_range, -rho))
    #     y_index_end = np.where(self.y_range == self._find_nearest(self.y_range, rho))
    #
    #     y_value = self.y_range[:int(y_index_end[0]) + 1]
    #     # X, Y = np.meshgrid(self.x_range, self.y_range)
    #
    #     b = np.ndarray([[self.get_b_abs((x, y, 0)) for x in x_value] for y in y_value])
    #     plt.imshow(b,  aspect='auto', extent=[start, end, -rho, rho])
    #     plt.colorbar()
    #     plt.show()

    def plot_2d_vectormap(self, start, end, rho=0):
        """Plot a 2D vectormap.

        WORKING FUNCTION
        """

        self.calculate_b_field(zero=start, meshsize=(end - start, rho, 0))

        x_index_start = np.where(self.x_range == self._find_nearest(self.x_range, start))
        x_index_end = np.where(self.x_range == self._find_nearest(self.x_range, end))

        x_value = self.x_range[int(x_index_start[0]):int(x_index_end[0]) + 1]
        # print(self.x_range)

        # y_index_start = np.where(self.y_range == self.find_nearest(self.y_range, -rho))
        y_index_end = np.where(self.y_range == self._find_nearest(self.y_range, rho))

        y_value = self.y_range[:int(y_index_end[0]) + 1]
        distance_of_arrows = int(len(y_value) / 25.)+1
        y_pos = y_value[::distance_of_arrows]
        x_pos = x_value[0::int(1. * len(x_value) / len(y_pos)) + 1]

        # print(x_pos, y_pos)

        b_vec = np.array([[self.b[(x, y, 0)] for x in x_value] for y in y_value])

        # print(B_vec)
        u = b_vec[:, :, 0]
        v = b_vec[:, :, 1]
        # print(x_value, y_value)

        b = np.array([[self.get_b_abs((x, y, 0)) for x in x_value] for y in y_value])
        # print(b_vec)
        # print(b)

        plot_range_min = - rho
        plot_range_max = + rho

        if not (plot_range_max - plot_range_min):
            plot_range_max += 0.1
            plot_range_min -= 0.1

        plt.imshow(b, aspect='auto', extent=[start, end, plot_range_min, plot_range_max])
        plt.colorbar()
        plt.quiver(x_pos, y_pos, u/np.sqrt(u**2+v**2), v/np.sqrt(u**2+v**2))
        plt.show()

    def get_b(self, arg):
        x, y, z = arg
        b = np.zeros(3)
        local_b = self.b[(x, abs(y), abs(z))]
        b[0] = local_b[0]
        b[1] = np.sign(y)*local_b[1]
        b[2] = np.sign(z)*local_b[2]

        return b

    def get_b_abs(self, arg):
        x, y, z = arg
        # b = np.zeros(3)
        local_b = self.b[(x, y, z)]
        return np.linalg.norm(local_b)
