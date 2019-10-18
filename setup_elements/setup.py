# -*- coding: utf-8 -*-

from setup_elements.elements import Coil, RealCoil, SquareCoil
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
import itertools
from scipy.interpolate import interp1d


class Setup:

    def __init__(self, increment, speed='fast'):
        self.elements = []
        self.B = None
        self.setup_changed = False
        self.x_range = None
        self.y_range = None
        self.z_range = None
        self.I = None
        self.increment = increment
        self.speed = speed

    def _find_nearest(self, array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return array[idx]

    def create_coil(self, coil_mid_pos=0, length=0.1, windings=1000, current=10,  r=0.05, wire_d=0.006, angle_y=0, angle_z = 0):
        if self.speed == 'fast':
            element = SquareCoil()
        else:
            element = RealCoil()

        element.create_coil(coil_mid_pos, length, windings, current,  r, wire_d, angle_y, angle_z)
        self.elements.append(element)
        self.setup_changed = True

    def b_x(self, x, rho=0):
        field = 0
        for element in self.elements:
            field += element.b_x(x, rho)
        return field

    def b_rho(self, x, rho):
        field = 0
        for element in self.elements:
            np.add(field, element.b_rho(x, rho), out=field)

        return field

    def B_field(self, r: '(x, y, z)'):
        field = np.array((0., 0., 0.))
        for element in self.elements:
            np.add(field, element.B_field(*r), out=field)
            # field += element.B_field(*r)

        return field

    def change_current(self, current):
        self.I = current
        self.setup_changed = True

    def calculate_b_field(self, zero, meshsize: '(x, y, z)'):
        
        if self.B is None or self.setup_changed:

            self.x_range = np.arange(zero, meshsize[0] + zero, self.increment)
            self.y_range = np.arange(-meshsize[1], meshsize[1]+self.increment, self.increment)
            self.z_range = np.arange(-meshsize[2], meshsize[2]+self.increment, self.increment)
            # print(self.x_range, self.y_range, self.z_range)
            args = list(itertools.product(self.x_range, self.y_range, self.z_range))
            print(len(args), " calculations")
            print('calculate')

            with Pool(4) as p:
                result = p.map(self.B_field, args)

            print('calculation finished')

            self.B = dict(zip(args, result))

        if meshsize[0]-1 > max(self.x_range) or meshsize[1]-1 > max(self.y_range) or  meshsize[2]-1 > max(self.z_range):

            self.x_range = np.arange(zero, meshsize[0] + zero, self.increment)
            self.y_range = np.arange(-meshsize[1], meshsize[1], self.increment)
            self.z_range = np.arange(-meshsize[2], meshsize[2], self.increment)
            args = list(itertools.product(self.x_range, self.y_range, self.z_range))

            for key in self.B.keys():
                if key in args:
                    args.remove(key)

            print('calculate extension')

            with Pool(4) as p:
                result = p.map(self.B_field, args)

            print('calculation finished')

            self.B.update(dict(zip(args, result)))

        self.setup_changed = False

    def return_1D_FI(self, start, end, rho=0):
        self.calculate_b_field(zero=start, meshsize=(end - start, rho, 0))

        x_index_start = np.where(self.x_range == self._find_nearest(self.x_range, start))
        x_index_end = np.where(self.x_range == self._find_nearest(self.x_range, end))

        x_values = self.x_range[int(x_index_start[0]):int(x_index_end[0])+1]

        if rho != 0:
            y_values = [np.linalg.norm(self.B[(x, rho, 0)]) for x in x_values]
        else:
            y_values = [self.B[(x, 0, 0)][0] for x in x_values]

        return sum(y_values)*self.increment




    def plot_1D_abs(self, start, end, rho=0):

        self.calculate_b_field(zero=start, meshsize=(end - start, rho, 0))

        x_index_start = np.where(self.x_range == self._find_nearest(self.x_range, start))
        x_index_end = np.where(self.x_range == self._find_nearest(self.x_range, end))

        x_values = self.x_range[int(x_index_start[0]):int(x_index_end[0]) + 1]

        if rho != 0:
            y_values = [np.linalg.norm(self.get_B((x, rho, 0))) for x in x_values]
        else:
            y_values = [self.B[(x, 0, 0)][0] for x in x_values]
        plt.plot(self.x_range, y_values)
        plt.show()

    def plot_1D_vector(self, start, end, rho=0):

        self.calculate_b_field(zero=start, meshsize=(end - start, rho, 0))
        x_index_start = np.where(self.x_range == self._find_nearest(self.x_range, start))
        x_index_end = np.where(self.x_range == self._find_nearest(self.x_range, end))

        x_pos = self.x_range[int(x_index_start[0]):int(x_index_end[0]) + 1]

        y_pos = 0
        B_vec = np.array([self.B[(x, rho, 0)][:2] for x in x_pos])
        u = B_vec[:,0]
        v = B_vec[:,1]

        plt.quiver(x_pos, y_pos, u/np.sqrt(u**2+v**2), v/np.sqrt(u**2+v**2))
        plt.show()

    def data_compare_1D(self, x_data, y_data, rho=0):

        self.calculate_b_field(zero=min(x_data), meshsize=(max(x_data), rho, 0))

        x_values = self.x_range
        #print(x_data, x_values)

        if rho != 0:
            y_values = [np.linalg.norm(self.B[(x, rho, 0)]) for x in x_values]

        else :
            y_values = [self.B[(x, rho, 0)][0] for x in x_values]

        B = interp1d(x_data, y_data)

        B_data = np.array([B(x) for x in x_values])
        plt.plot(x_values, y_values-B_data)
        plt.show()

    def plot_2D_vector(self, start, end, rho=0):

        self.calculate_b_field(zero=start, meshsize=(end - start, rho, 0))
        x_index_start = np.where(self.x_range == self._find_nearest(self.x_range, start))
        x_index_end = np.where(self.x_range == self._find_nearest(self.x_range, end))

        x_value = self.x_range[int(x_index_start[0]):int(x_index_end[0]) + 1]
        distance_of_arrows = int(len(self.y_range)/25.)+1

        y_pos = self.y_range[::distance_of_arrows]

        x_pos = x_value[0::int(1.*len(x_value)/len(y_pos))+1]

        B_vec = np.array([[self.B[(x,y,0)] for x in x_pos] for y in y_pos])
        #print(B_vec)
        u = B_vec[:,:,0]
        v = B_vec[:,:,1]

        plt.quiver(x_pos, y_pos, u/np.sqrt(u**2+v**2), v/np.sqrt(u**2+v**2))
        plt.show()

    def plot_2D_map(self, start, end, rho):

        self.calculate_b_field(zero=start, meshsize=(end - start, rho, 0))
        x_index_start = np.where(self.x_range == self._find_nearest(self.x_range, start))
        x_index_end = np.where(self.x_range == self._find_nearest(self.x_range, end))

        x_value = self.x_range[int(x_index_start[0]):int(x_index_end[0]) + 1]

        #y_index_start = np.where(self.y_range == self.find_nearest(self.y_range, -rho))
        y_index_end = np.where(self.y_range == self._find_nearest(self.y_range, rho))

        y_value = self.y_range[:int(y_index_end[0]) + 1]
        #X, Y = np.meshgrid(self.x_range, self.y_range)
        B = np.matrix([[self.get_B_abs((x, y, 0)) for x in x_value] for y in y_value])
        plt.imshow(B,  aspect='auto', extent=[start, end, -rho, rho])
        plt.colorbar()
        plt.show()

    def plot_2D_vectormap(self, start, end, rho=0):

        self.calculate_b_field(zero=start, meshsize=(end - start, rho, 0))

        x_index_start = np.where(self.x_range == self._find_nearest(self.x_range, start))
        x_index_end = np.where(self.x_range == self._find_nearest(self.x_range, end))

        x_value = self.x_range[int(x_index_start[0]):int(x_index_end[0]) + 1]
        #print(self.x_range)

        #y_index_start = np.where(self.y_range == self.find_nearest(self.y_range, -rho))
        y_index_end = np.where(self.y_range == self._find_nearest(self.y_range, rho))

        y_value = self.y_range[:int(y_index_end[0]) + 1]
        # print(x_pos)
        distance_of_arrows = int(len(y_value) / 25.)+1
        y_pos = y_value[::distance_of_arrows]
        x_pos = x_value[0::int(1. * len(x_value) / len(y_pos)) + 1]
        B_vec = np.array([[self.B[(x,y,0)] for x in x_pos] for y in y_pos])
        # print(B_vec)
        u = B_vec[:, :, 0]
        v = B_vec[:, :, 1]
        B = np.matrix([[self.get_B_abs((x, y, 0)) for x in x_value] for y in y_value])
        plt.imshow(B, aspect='auto', extent=[start, end, -rho, rho])
        plt.colorbar()
        plt.quiver(x_pos, y_pos, u/np.sqrt(u**2+v**2), v/np.sqrt(u**2+v**2))
        plt.show()

    # def get_B(self, arg):
    #     x, y, z = arg
    #     B = np.zeros(3)
    #     Local_B = self.B[(x, abs(y), abs(z))]
    #     B[0] = Local_B[0]
    #     B[1] = np.sign(y)*Local_B[1]
    #     B[2] = np.sign(z)*Local_B[2]
    #
    #     return B

    def get_B_abs(self, arg):
        x, y, z = arg
        # B = np.zeros(3)
        local_b = self.B[(x, y, z)]
        B = local_b
        return np.linalg.norm(B)