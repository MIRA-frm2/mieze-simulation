# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""General setup allowing placement of different elements."""

import itertools
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from multiprocessing import Pool
import numpy as np

from elements.coils import Coil, RealCoil, SquareCoil
from utils.helper_functions import read_data_from_file


class Setup:
    """Class that simulates a physical setup."""
    def __init__(self, increment):
        self.elements = []
        self.b = None
        self.b_cartesian = None

        self.setup_changed = False

        # Computational discretized space
        self.x_range = None
        self.y_range = None
        self.z_range = None

        self.start = None
        self.end = None
        self.rho = None

        self.current = None
        self.increment = increment

    def create_setup(self, current, l1=None, l2=None):
        raise NotImplementedError

    def create_element(self, coil_type, coil_mid_pos=0, length=0.1, windings=1000, current=10,  r=0.05,
                       wire_d=0.006, angle_y=0, angle_z=0):
        """Create the physical geometry of the coils."""
        self.elements.append(coil_type(coil_mid_pos, length, windings, current,  r, wire_d, angle_y, angle_z))
        self.setup_changed = True

    @staticmethod
    def _find_nearest(array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value).argmin())
        return idx

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
        """Change the current value."""
        self.current = current
        self.setup_changed = True

    def _get_square_element(self):
        return self.elements[0]

    def calculate_b_field(self, coordinate_system='cartesian', **kwargs):
        """Calculate the magnetic field."""

        if coordinate_system == 'cartesian':
            square_element = self._get_square_element()

            zoom_factor = 2

            start = kwargs.pop('start', -0.5)
            plane = kwargs.pop('plane', True)
            # coil_distance = kwargs.pop('coil_distance', square_element.width/10)

            if self.b is None or self.setup_changed:
                if plane == 'yz':
                    self.x_range = np.arange(1)
                else:
                    self.x_range = np.arange(start, square_element.length + self.increment, self.increment)

                if plane == 'xz':
                    self.y_range = np.arange(1)
                else:
                    self.y_range = np.arange(- zoom_factor * square_element.width - self.increment/2,
                                             zoom_factor * square_element.width + self.increment/2,
                                             self.increment)

                if plane == 'xy':
                    self.z_range = np.arange(1)
                else:
                    self.z_range = np.arange(- zoom_factor * square_element.height - self.increment/2,
                                             zoom_factor * square_element.height + self.increment/2,
                                             self.increment)

                # print(self.x_range, self.y_range, self.z_range)
                args = list(itertools.product(self.x_range, self.y_range, self.z_range))
                print(len(args), " calculations")
                print('calculate')

                with Pool(4) as p:
                    result = p.map(self.b_field, args)

                print('calculation finished')

                self.b = dict(zip(args, result))
                # print(self.b)

        elif coordinate_system == 'cylindrical':
            start = kwargs.pop('start', 0)
            end = kwargs.pop('end', 0)
            rho = kwargs.pop('rho', 0)
            plane = kwargs.pop('plane', True)

            resolution_factor = 1.125
            compute_range = resolution_factor * rho

            if plane == 'yz':
                self.x_range = np.arange(1)
            else:
                self.x_range = np.arange(start, end - start, self.increment)

            if plane == 'xz':
                self.y_range = np.arange(1)
            else:
                self.y_range = np.arange(- compute_range, compute_range, self.increment)

            if plane == 'xy':
                self.z_range = np.arange(1)
            else:
                self.z_range = np.arange(- compute_range, compute_range, self.increment)

            args = list(itertools.product(self.x_range, self.y_range, self.z_range))

            if self.b is None or self.setup_changed:
                print(len(args), " calculations")
                print('calculate')

                with Pool(4) as p:
                    result = p.map(self.b_field, args)

                print('calculation finished')

                self.b = dict(zip(args, result))
                print(self.b)

            # if meshsize[0]-1 > max(self.x_range) or meshsize[1]-1 > max(self.y_range)
            # or meshsize[2]-1 > max(self.z_range):
            #
            #     self.x_range = np.arange(zero, meshsize[0] + zero, self.increment)
            #     self.y_range = np.arange(-meshsize[1], meshsize[1], self.increment)
            #     self.z_range = np.arange(-meshsize[2], meshsize[2], self.increment)
            #     args = list(itertools.product(self.x_range, self.y_range, self.z_range))
            #
            #     for key in self.b.keys():
            #         if key in args:
            #             args.remove(key)
            #
            #     print('calculate extension')
            #
            #     with Pool(4) as p:
            #         result = p.map(self.b_field, args)
            #
            #     print('calculation finished')
            #
            #     self.b.update(dict(zip(args, result)))

        self.setup_changed = False

    def get_plot_data(self, plane, source='new'):
        if source == 'storage':
            b, extent = read_data_from_file()
        else:
            b = self.get_magnetic_field_value(plane, plane_position=0)

            if plane == 'yz':
                extent = (self.y_range.min(), self.y_range.max(), self.z_range.min(), self.z_range.max())
            elif plane == 'xy':
                extent = (self.x_range.min(), self.x_range.max(), self.y_range.min(), self.y_range.max())
            elif plane == 'xz':
                extent = (self.x_range.min(), self.x_range.max(), self.z_range.min(), self.z_range.max())
            else:
                extent = None
        return b, extent

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

    def plot_field_1d_abs(self):

        # if self.rho != 0:
        #     y_values = [[[np.linalg.norm(self.get_b_abs((x, y, z))) for x in self.x_range] for y in self.y_range]
        #     for z in self.z_range][0][0]
        # else:
        #     try:
        #         y_values = [self.b[(x, 0, 0)][0] for x in self.x_range]
        #     except KeyError:
        #         # Todo: How to plot b field for non-zero y and z values?
        #         y_values = None

        y_values = [[[np.linalg.norm(self.get_b_abs((x, y, z))) for x in self.x_range] for y in self.y_range] for z in
                    self.z_range][0][0]

        print(self.x_range)
        print(y_values)

        plt.plot(self.x_range, y_values)
        plt.show()

    def plot_field_2d_abs(self, plane):
        b, extent = self.get_plot_data(plane)

        plt.imshow(b, aspect='auto', cmap=cm.magma, extent=extent)

        plt.colorbar()
        plt.show()

    def get_magnetic_field_value(self, plane, plane_position):
        """Return the magnetic field value at a given plane."""

        if plane == 'xy':
            component = 2
            idx = self._find_nearest(self.z_range, plane_position)
        elif plane == 'xz':
            component = 1
            idx = self._find_nearest(self.y_range, plane_position)
        elif plane == 'yz':
            component = 0
            idx = self._find_nearest(self.x_range, plane_position)
        else:
            component = None
            idx = None

        b = np.array([[[self.b[(x, y, z)][component] for z in self.z_range]
                       for y in self.y_range] for x in self.x_range])

        if component == 0:
            return b[idx]
        elif component == 1:
            return b[:, idx]
        if component == 2:
            return b[:, :, idx]

    def get_b(self, arg):
        # ToDo: why absolute values of y and z
        # Legacy:
        # x, y, z = arg
        # b = np.zeros(3)
        # print(self.b)
        # local_b = self.b[(x, abs(y), abs(z))]
        # b[0] = local_b[0]
        # b[1] = np.sign(y)*local_b[1]
        # b[2] = np.sign(z)*local_b[2]
        # return b
        x, y, z = arg
        local_b = self.b[(x, y, z)]
        return local_b

    def get_b_abs(self, arg):
        """

        Parameters
        ----------
        arg: tuple
            x, y, z coordinates
        """
        # b = np.zeros(3)
        # print(self.b)
        local_b = self.b[arg]
        return np.linalg.norm(local_b)
