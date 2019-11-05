# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""An implementation of the polariser."""

from experiments.mieze.parameters import polarizer_position
import numpy as np
# from scipy.optimize import curve_fit


class polariser:
    MU_0 = 4e-7 * np.pi  # N/A, vacuum permeability

    def __init__(self):
        pass
    c = 0.05  # shift polarisator position
    m = np.array((0, 90, 0))  # magnetic dipol moment Polarisator; both obtained by real measurements and data analysis


    @classmethod
    def B_field(cls, x, y, z, zero=polarizer_position):
        x -= zero
        r_vec = np.array((x+cls.c, y, z))
        r = np.linalg.norm(r_vec)
        B = (cls.MU_0/(4*np.pi*r**2))*((3*r_vec*np.linalg.norm(cls.m*r_vec) - cls.m*r**2) / r**3)
        # print()
        # B = cls.MU_0/(4*np.pi*(r+cls.c)**3)*cls.m
        # B = (1e-7)*cls.m/(x+cls.c)**3
        return B*1e4
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
