# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""General setup allowing placement of different elements."""

import itertools
import logging
from multiprocessing import Pool
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm


from utils.helper_functions import read_data_from_file, find_nearest, add_earth_magnetic_field

# Create a custom logger
logger = logging.getLogger(__name__)


class Setup:
    """Class that simulates a physical setup."""
    def __init__(self, consider_earth_field=False):
        self.elements = []
        self.b = None
        self.b_cartesian = None

        self.setup_changed = False

        # Computational discretized space
        self.x_range = None
        self.y_range = None
        self.z_range = None

        self.mesh = None

        self.start = None
        self.end = None
        self.rho = None

        self.current = None

        self.x_ticks = list()
        self.x_ticks_labels = list()

        self.computation_number = 0

        self.consider_earth_field = consider_earth_field

    def create_setup(self):
        raise NotImplementedError

    def create_element(self, element_class, position, **kwargs):
        """Create the physical geometry of the coils."""
        self.elements.append(
            element_class(position=position, **kwargs))

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
        self.computation_number += 1

        field = np.array((0., 0., 0.))

        for element in self.elements:
            np.add(field, element.b_field(*r), out=field)

        # logger.error(f'{self.computation_number}: Computed total magnetic field for point {r} is {field}')

        field = add_earth_magnetic_field(field, flag=self.consider_earth_field)

        return field

    def change_current(self, current):
        """Change the current value."""
        self.current = current
        self.setup_changed = True

    def _get_square_element(self):
        return self.elements[0]

    def initialize_computational_space(self, **kwargs):
        """Initialize the 3d discretized computational space."""
        x_start = kwargs.pop('x_start', 0)
        x_end = kwargs.pop('x_end', 1)
        x_step = kwargs.pop('x_step', 0.1)

        y_start = kwargs.pop('y_start', 0)
        y_end = kwargs.pop('y_end', 1)
        z_start = kwargs.pop('z_start', 0)
        z_end = kwargs.pop('z_end', 1)
        yz_step = kwargs.pop('yz_step', 0.1)

        self.x_range = np.arange(x_start, x_end + x_step, x_step)
        self.y_range = np.arange(y_start, y_end + yz_step, yz_step)
        self.z_range = np.arange(z_start, z_end + yz_step, yz_step)

    def calculate_b_field(self, point=None):
        """Calculate the magnetic field."""
        if point:
            return self.b_field(point)
        else:
            # logger.error(f'{self.x_range}, {self.y_range}, {self.z_range}')
            positions = list(itertools.product(self.x_range, self.y_range, self.z_range))
            logger.error(f'{len(positions)} calculations')
            logger.info('calculate')

            with Pool(4) as p:
                result = p.map(self.b_field, positions)

            logger.info('calculation finished')

            self.b = dict(zip(positions, result))

    def get_b_vec(self):
        """Returns the magnetic field as a vector."""
        bx = [self.b[x, 0, 0][0] for x in self.x_range]
        by = [self.b[x, 0, 0][1] for x in self.x_range]
        bz = [self.b[x, 0, 0][2] for x in self.x_range]
        return bx, by, bz

    def _get_plane_position(self, component, plane_position):
        """Get the plane position from the 3D grid."""
        if component == 'z':
            component = 2
            plane_idx = find_nearest(self.z_range, plane_position)
        elif component == 'y':
            component = 1
            plane_idx = find_nearest(self.y_range, plane_position)
        elif component == 'x':
            component = 0
            plane_idx = find_nearest(self.x_range, plane_position)
        else:
            component = 'abs'
            plane_idx = find_nearest(self.y_range, plane_position)
        return component, plane_idx

    @staticmethod
    def _get_b_field_values_from_plane(b, component, plane_idx):
        if component == 0:
            return b[plane_idx]
        elif component == 1:
            return b[:, plane_idx]
        if component == 2:
            return b[:, :, plane_idx]

    def _get_b_field_values(self):
        return np.array([[[self.b[(x, y, z)] for z in self.z_range]
                        for y in self.y_range] for x in self.x_range])

    def get_magnetic_field_value(self, component, plane_position):
        """Return the magnetic field value at a given plane."""

        component, plane_idx = self._get_plane_position(component, plane_position)

        if component == 'abs':
            b = self._get_b_field_values()
        else:
            b = self._get_b_field_values()[component]

        return self._get_b_field_values_from_plane(b, component, plane_idx)

    def get_2d_abs_plot_data(self, plane, source='new'):
        if source == 'storage':
            b, extent = read_data_from_file('../data/data.csv')
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

    def get_1d_b_values(self, component):
        return self.get_magnetic_field_value(component='abs', plane_position=0)

    @staticmethod
    def get_numerical_component(component):
        if component == 'x':
            return 0
        elif component == 'y':
            return 1
        elif component == 'z':
            return 2

    def plot_field_1d_scalar(self, component):
        component = self.get_numerical_component(component)

        plot_x_values = list()
        plot_y_values = list()

        for point, b_field in self.b.items():
            x_value = point[0]
            if x_value not in plot_x_values:
                plot_x_values.append(x_value)

            y_value = b_field[component]
            if y_value not in plot_y_values:
                plot_y_values.append(y_value)

        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])  # main axes

        ax.plot(plot_x_values, plot_y_values)

        locs, labels = plt.xticks()
        if self.x_ticks and self.x_ticks_labels:
            print(f'{locs}, {labels}, {self.x_ticks}, {self.x_ticks_labels}')
            plt.xticks([*locs, *self.x_ticks], [*labels, *self.x_ticks_labels])
            print(f'{plt.xticks()}')

        plt.xlabel('Beam line [m]')
        plt.ylabel('Magnetic field [G]')

        plt.show()

    def plot_field_1d_vec(self, _type='3d'):
        # Redefine plot ranges
        # y = len(self.x_range) * [0]
        # z = len(self.x_range) * [0]
        #
        # bx, by, bz = self.get_b_vec()
        # print(f'bx:{bx}\nby:{by}\nbz_{bz}')
        # x, y, z = np.meshgrid(self.x_range, self.y_range, self.z_range)

        if _type == '3d':
            fig = plt.figure()

            ax = fig.add_subplot(111, projection='3d')

            for point, b_field in self.b.items():
                ax.quiver(point[0], point[1], point[2],
                          b_field[0], b_field[1], b_field[2],
                          # pivot='tip', length=vlength, arrow_length_ratio=0.3/vlength
                          )

            ax.set_xlim([min(self.x_range), max(self.x_range)/3])
            ax.set_ylim([-0.01, 0.01])
            ax.set_zlim([-0.0025, 0.0025])

            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_zlabel('z')

        elif _type == '2d':
            xi, yi = np.meshgrid(self.x_range, 0, indexing='ij')
            for x in self.x_range:
                v = self.b[x, 0, 0]
                print(v)
                # plt.axes([0.065, 0.065, 0.9, 0.9])
                plt.quiver(xi, yi, v[0], v[2], alpha=.5)
                plt.quiver(xi, yi, v[0], v[2], edgecolor='k', facecolor='none', linewidth=.5)

        plt.show()

    def plot_field_2d_abs(self, plane):
        b, extent = self.get_2d_abs_plot_data(plane)

        plt.imshow(b, aspect='auto', cmap=cm.magma, extent=extent)

        plt.colorbar()
        plt.show()

    def get_2d_vec_plot_data(self):
        by = np.array([[[self.b[(x, y, z)][0] for z in self.z_range]
                       for y in self.y_range] for x in self.x_range])[:, 0]

        bz = np.array([[[self.b[(x, y, z)][1] for z in self.z_range]
                       for y in self.y_range] for x in self.x_range])[:, :, 0]
        return by, bz

    def plot_field_2d_vec(self):
        by, bz = self.get_2d_vec_plot_data()

        plt.quiver(self.y_range, self.z_range, by, bz)
        plt.xlim(-1, 1)
        plt.ylim(-1, 1)

        plt.show()

    def set_plot_ticks(self, set_ticks=False):
        if set_ticks:
            for element in self.elements:
                self.x_ticks.append(element.position_x)
                self.x_ticks_labels.append(element.name)