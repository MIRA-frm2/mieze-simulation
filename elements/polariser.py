# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""An implementation of the polariser."""

import logging
import numpy as np
# from scipy.optimize import curve_fit2

from elements.base import BasicElement
from experiments.mieze.parameters import POLARISATOR

from utils.physics_constants import MU_0
from utils.helper_functions import get_vector_norm


# Create a custom logger
logger = logging.getLogger(__name__)


class Polariser(BasicElement):

    def __init__(self, position=(POLARISATOR, 0, 0), *args, **kwargs):

        super(Polariser, self).__init__(position, name='polariser')

        self.position = position

        self.c = kwargs.get('c', 0.05)  # shift polarisator position
        # magnetic dipol moment Polarisator; both obtained by real measurements and data analysis
        self.m = kwargs.get('m', np.array((0, 90, 0)))

    def b_field(self, x, y, z):
        x -= self.position_x
        r_vec = np.array((x + self.c, y, z))
        r = get_vector_norm(r_vec)

        prefactor = (MU_0 / (4 * np.pi * r**2))
        elem1 = 3 * r_vec * get_vector_norm(self.m * r_vec)
        elem2 = self.m * r**2

        b = prefactor * (elem1 - elem2) / r**3
        # logger.error(f'{b}')
        # B = cls.MU_0/(4*np.pi*(r+cls.c)**3)*cls.m
        # B = (1e-7)*cls.m/(x+cls.c)**3
        return b[0] * 1e4 * b
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
