# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""An implementation of the Spin Flipper."""

from experiments.mieze.parameters import SpinFlipper_position1, SpinFlipper_position2
from numpy import arctan

from elements.base import BasicElement

from utils.physics_constants import MU_0, pi, unit


class SpinFlipper(BasicElement):

    def __init__(self, position=None, *args, **kwargs):
        super(SpinFlipper, self).__init__()

        self.windings = kwargs.get('windings', 100)
        self.length = kwargs.get('length', 13e-2 * unit.m)  # [m], winding length of the coil
        self.thickness = kwargs.get('thickness', 1e-2 * unit.m)  # [m], thickness of the coil
        self.current = kwargs.get('current', 1.0 * unit.ampere)  # [A], coil currents

    def _sf_th(self, current, position):
        x0 = self.length / 2.0
        n = self.windings / self.length

        b_eff = n * MU_0 * current / pi * arctan(x0 / position) * unit.m / unit.radian
        return b_eff
    
    def b_field(self, x, y, z):
        zero = SpinFlipper_position1
        # if sf_name == 'sf1':
        #     zero = SpinFlipper_position1
        # elif sf_name == 'sf2':
        #     zero = SpinFlipper_position2
        # else:
        #     raise RuntimeError('Wrong name for sf given.')

        y = x - zero

        b1 = self._sf_th(self.current, y + self.thickness / 2.0)
        b2 = -self._sf_th(self.current, y - self.thickness / 2.0)
        return (b1 + b2) * 1e4
