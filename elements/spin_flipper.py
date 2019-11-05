# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""An implementation of the Spin Flipper."""

import numpy as np
# from experiments.mieze.parameters import SpinFlipper_position


class SpinFlipper:
    MU_0 = 4e-7 * np.pi  # N/A, vacuum permeability
    N = 100
    LENGTH = 13e-2  # m, winding length of the coil
    THICKNESS = 1e-2  # m, thickness of the coil

    def __init__(self, position):
        self.position = position

    @classmethod
    def _sf_th(cls, current, position):
        x0 = cls.LENGTH / 2.0
        n = cls.N / cls.LENGTH

        B_eff = n * cls.MU_0 * current / np.pi * np.arctan(x0 / position)
        return B_eff
    
    @classmethod
    def sf(cls, sf_name:str, current:float, z_position):
    
        if sf_name == 'sf1':
            zero = Positions.get_position_sf1()
        elif sf_name == 'sf2':
            zero = Positions.get_position_sf2()
        else:
            raise RuntimeError('Wrong name for sf given.')

        y = z_position - zero

        B1 = cls._sf_th(current, y + cls.THICKNESS / 2.0)
        B2 = -cls._sf_th(current, y - cls.THICKNESS / 2.0)
        return (B1 + B2) * 1e4



        

