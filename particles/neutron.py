# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Neutron class."""


import itertools
from multiprocessing import Pool
import numpy as np
import numpy.random as r
import os
import random


from utils.helper_functions import save_obj, load_obj

cwd = os.getcwd()


class Neutron:
    """Implements neutrons and its properties."""

    def __init__(self, speed, polarisation, position):
        self.speed = speed
        self.polarisation = polarisation
        self.position = position

    @staticmethod
    def _omega(b):
        gamma = 1.83247172e4
        return b*gamma

    @staticmethod
    def _precession_angle(time, b):
        gamma = 1.83247172e4
        return gamma * np.linalg.norm(b)*time

    @staticmethod
    def _rotate(vektor, phi, axis):
        n = axis/np.linalg.norm(axis)
        c = np.cos(phi)
        s = np.sin(phi)

        n1 = n[0]
        n2 = n[1]
        n3 = n[2]

        R = [[n1**2 * (1-c) + c, n1 * n2 * (1-c) - n3 * s, n1 * n3 * (1-c) + n2 * s],
             [n2 * n1 * (1-c) + n3 * s, n2**2 * (1-c) + c, n2 * n3 * (1-c) - n1 * s],
             [n3 * n1 * (1-c) - n2 * s, n3 * n2 * (1-c) + n1 *s, n3**2*(1-c) + c]]
        return np.dot(R, vektor)
    
    def _polarisation_change(self, neutron, b, L):
        t = self._time_in_field(velocity=neutron['speed'])
        phi = self._precession_angle(t, b)
        return self._rotate(vektor=neutron['polarisation'], phi=phi, axis=b)
    
    def create_b_map(self, b_function, profiledimension):
        """

        Parameters
        ----------
        b_function
        profiledimension: tuple

        Returns
        -------

        """
        if not len(profiledimension) == 2:
            print('error: profile dimension has to be 2D ')
            return -1
        
        max_y = profiledimension[0] 
        max_z = profiledimension[1]
        
        x = np.round(np.arange(0, self.totalflightlength, self.incrementsize),
                     decimals=-int(np.log10(self.incrementsize)))
        y = np.round(np.arange(-max_y, max_y + self.incrementsize, self.incrementsize),
                     decimals=-int(np.log10(self.incrementsize)))
        z = np.round(np.arange(-max_z, max_z + self.incrementsize, self.incrementsize),
                     decimals=-int(np.log10(self.incrementsize)))

        args = list(itertools.product(x, y, z))
        
        with Pool(4) as p:
            result = p.map(b_function, args)

        bmap = dict(zip(args, result))
        
        # self.B_map = bmap
        save_obj(bmap, 'B_map')
    
    def simulate_neutrons(self):
        toberemoved = []

        b_map = load_obj('../data/data')
        # try:
        #     pass
        # except FileNotFoundError:
        #     # Todo: If not magnetic field is found, create one.
        #     # self.create_b_map(b_function, (0.001, 0.001))
        #     # b_map = self.load_obj('B_map')
        #     pass

        for neutron in self.neutrons:
            # ToDo: refactor
            # check if it is in the calculated beamprofile (y,z plane)
            if not (0, neutron['beampositiony'], neutron['beampositionz']) in b_map:
                toberemoved.append(neutron)
                continue

            for ii in range(int(self.totalflightlength/self.incrementsize)):
                i = round(ii*self.incrementsize, ndigits=3)
                neutron['polarisation'] = self._polarisation_change(neutron, b_map[(i, neutron['beampositiony'], neutron['beampositionz'])], self.incrementsize)

        for neutron in toberemoved:
            self.neutrons.remove(neutron)
                
    def get_pol(self):
        return self.polarisation

    def reset_pol(self):
        c = 0.31225  # sqrt(1-polarisierung²)
        x = c*r.rand()
        z = np.sqrt(c**2-x**2)

        self.polarisation = np.array([x, 0.95, z])
