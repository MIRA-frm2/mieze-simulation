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
import pickle
import random


cwd = os.getcwd()


class Neutrons:
    """Implements neutrons and its properties."""

    def __init__(self, number_of_neutrons, velocity = 3956/4.5, incrementsize=0.001, totalflightlength=3, beamsize=0.02):
        self.neutrons = []
        self.number_of_neutrons = number_of_neutrons
        self.velocity = velocity
        self.incrementsize = incrementsize
        self.totalflightlength = totalflightlength
        while len(self.neutrons)<self.number_of_neutrons:
            c = 0.31225 # sqrt(1-polarisierung²)
            x = c*r.rand()
            z = np.sqrt(c**2-x**2)
            speed= velocity
            neutron = {}
            neutron['speed'] = random.gauss(speed, 0.02*speed)
            neutron['polarisation'] = np.array([x, 0.95, z])
            neutron['beampositiony'] = round(random.gauss(0, beamsize/5), ndigits=-int(np.log10(self.incrementsize)))
            neutron['beampositionz'] = round(random.gauss(0, beamsize/5), ndigits=-int(np.log10(self.incrementsize)))
        
            self.neutrons.append(neutron)

    @staticmethod
    def save_obj(obj, name):
        with open(f'obj/{name}.pkl', 'wb') as f:
            pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

    @staticmethod
    def load_obj(name):
        with open(f'obj/{name}.pkl', 'rb') as f:
            return pickle.load(f)

    def _timeInField(self, velocity):
        return self.incrementsize/velocity
    
    def _omega(B):
        gamma = 1.83247172e4
        return B*gamma
    
    def _precessionAngle(self, time, B):
        gamma = 1.83247172e4
        return gamma*np.linalg.norm(B)*time
    
    def _rotate(self, vektor, phi, axis):
        n = axis/np.linalg.norm(axis)
        c = np.cos(phi)
        s = np.sin(phi)
        n1 = n[0]
        n2 = n[1]
        n3 = n[2]
        R = [[n1**2 *(1-c) + c , n1 *n2 *(1-c) - n3 *s , n1 *n3 *(1-c) + n2 *s ],[n2 *n1 *(1-c) + n3 *s , n2**2*(1-c) + c , n2 *n3 *(1-c) - n1 *s ],[n3 *n1 *(1-c) - n2 *s , n3 *n2 *(1-c) + n1 *s , n3**2*(1-c) + c]]
        return np.dot(R, vektor)
    
    def _polarisationChange(self, neutron, B, L):
        t = self._timeInField(velocity=neutron['speed'])
        phi = self._precessionAngle(t, B) 
        return self._rotate(vektor=neutron['polarisation'], phi=phi, axis=B)
    
    def create_B_map(self, B_function, profiledimension:tuple):
        
        if not len(profiledimension) == 2 :
            print('error: profile dimension has to be 2D ')
            return -1
        
        max_y = profiledimension[0] 
        max_z = profiledimension[1]
        
        x = np.round(np.arange(0, self.totalflightlength, self.incrementsize), decimals=-int(np.log10(self.incrementsize)))
        y = np.round(np.arange(-max_y, max_y+ self.incrementsize, self.incrementsize), decimals=-int(np.log10(self.incrementsize)))
        z = np.round(np.arange(-max_z, max_z+ self.incrementsize, self.incrementsize), decimals=-int(np.log10(self.incrementsize)))

        args = list(itertools.product(x,y,z))
        
        Bmap = {}

        with Pool(4) as p:
            result = p.map(B_function, args)

        Bmap = dict(zip(args, result))
        
        #self.B_map = Bmap
        self.save_obj(Bmap, 'B_map')
    
    def simulate_neutrons(self):
        toberemoved = []

        try:
            Bmap = self.load_obj('B_map')
        except FileNotFoundError:

            self.create_B_map(b_function, (0.001, 0.001))

        for neutron in self.neutrons:
            #check if it is in the calculated beamprofile (y,z plane)
            if not (0,neutron['beampositiony'], neutron['beampositionz']) in Bmap:
                toberemoved.append(neutron)
                continue

            for ii in range(round(self.totalflightlength/self.incrementsize)):
                i = round(ii*self.incrementsize, ndigits=3)
                neutron['polarisation'] = self._polarisationChange(neutron, Bmap[(i,neutron['beampositiony'], neutron['beampositionz'])], self.incrementsize)


        for neutron in toberemoved:
            self.neutrons.remove(neutron)
                
    def get_pol(self):
        
        xpol = np.average([neutron['polarisation'][0] for neutron in self.neutrons ])
        
        ypol = np.average([neutron['polarisation'][1] for neutron in self.neutrons ])
        
        zpol = np.average([neutron['polarisation'][2] for neutron in self.neutrons ])
        
        return xpol, ypol, zpol
        
    def reset_pol(self):
            c = 0.31225 # sqrt(1-polarisierung²)
            x = c*r.rand()
            z = np.sqrt(c**2-x**2)
            for neutron in self.neutrons:
                neutron['polarisation'] = np.array([x, 0.95, z])
