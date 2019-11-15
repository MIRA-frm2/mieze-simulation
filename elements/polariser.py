# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""An implementation of the polariser."""

import numpy as np
# from scipy.optimize import curve_fit2

from elements.base import BasicElement
from experiments.mieze.parameters import polarizer_position, unit

from utils.physics_constants import MU_0


class Polariser(BasicElement):

    def __init__(self, position=polarizer_position, *args, **kwargs):

        super(Polariser, self).__init__()

        self.position = position

        self.c = kwargs.get('c', 0.05 * unit.m)  # shift polarisator position
        # magnetic dipol moment Polarisator; both obtained by real measurements and data analysis
        self.m = kwargs.get('m', np.array((0 * unit.m, 90 * unit.m, 0 * unit.m)))

    def b_field(self, x, y, z):
        x *= unit.m
        y *= unit.m
        z *= unit.m

        x -= self.position
        r_vec = np.array((x + self.c, y, z))
        r = np.linalg.norm(r_vec)

        prefactor = (MU_0 / (4 * np.pi * r**2))
        elem1 = 3 * r_vec * np.linalg.norm(self.m * r_vec)
        elem2 = self.m * r**2
        b = prefactor * (elem1 - elem2) / r**3

        # print()
        # B = cls.MU_0/(4*np.pi*(r+cls.c)**3)*cls.m
        # B = (1e-7)*cls.m/(x+cls.c)**3
        return b*1e4
        # return 1e3*(cls.MU_0/4*np.pi)*cls.m / (x+cls.c)**3


# =============================================================================
#         def _mfd(x):
#             length = [0, 0.025, 0.05, 0.075, 0.10, 0.20, 0.30] #m
#             B = [52.5, 24, 9, 3.7, 2.1, 0.4, 0.12] #mT
#             length = np.array(length)
#             B = 10.0 * np. array(B)
#             popt, pcov = curve_fit(_polariser_mag, length, B, p0=(-5,1,0.1))
#             return _polariser_mag(x, *popt)
# 
# =============================================================================
